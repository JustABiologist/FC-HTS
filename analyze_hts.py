import argparse
import os
import re
import sys
import glob
import warnings

import flowio
import pandas as pd
import numpy as np

# Use non-interactive backend for compatibility (Windows, headless servers)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from openpyxl import load_workbook

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

def parse_arguments():
    parser = argparse.ArgumentParser(description="BD Fortessa HTS Analysis Pipeline")
    parser.add_argument("layout_file", help="Path to the Excel file containing the 96-well plate layout (names or BLANK).")
    parser.add_argument("fcs_input", help="Path to the directory containing .fcs files (or a single .fcs file).")
    parser.add_argument("wt_name", help="Name of the Wild Type (WT) annotation in the layout.")
    parser.add_argument("blank_name", help="Name of the Blank annotation in the layout used for background correction.")
    parser.add_argument("--channel", default="Blue C-A", help="Name of the channel to analyze (default: Blue-CA).")
    parser.add_argument("--flow-rate", type=float, default=None, help="Flow rate (µL/sec) for the instrument.")
    parser.add_argument("--doublet-threshold", type=float, default=1.5,
                        help="FSC-A / FSC-H ratio above which an event is treated as a doublet (default: 1.5).")
    parser.add_argument("--od-calibration", type=float, default=8e8,
                        help="Cells/mL per OD600 unit for your organism (default: 8e8 for E. coli in LB).")
    parser.add_argument("--cell-volume", type=float, default=20.0,
                        help="Volume of cell suspension added to each well in µL (default: 20).")
    parser.add_argument("--well-volume", type=float, default=300.0,
                        help="Total volume in each well in µL (default: 300).")
    parser.add_argument("--output", "-o", default=".", help="Output directory for results (default: current directory).")
    return parser.parse_args()

def clean_name(txt: str) -> str:
    """
    Remove leading/trailing whitespace, convert NB-space (chr(160)) to
    normal space and collapse runs of spaces to a single space.
    """
    if txt is None:
        return ''
    # replace NB-space by normal space
    txt = str(txt).replace('\xa0', ' ')
    # strip, collapse multiple spaces, make a single canonical version
    txt = re.sub(r'\s+', ' ', txt).strip()
    return txt

def normalize_well_id(well_str):
    """Normalizes well ID to A01 format."""
    match = re.match(r"([A-H])(\d+)", well_str.upper())
    if match:
        row = match.group(1)
        col = int(match.group(2))
        return f"{row}{col:02d}"
    return None

def read_plate_layout(excel_path):
    """
    Reads a 96-well plate layout from an Excel file.
    Supports 8x12 or partial grids (e.g. 6x12) without headers.
    Returns a dictionary mapping 'A01' -> 'SampleName'.
    """
    print(f"Reading plate layout from {excel_path}...")
    try:
        df = pd.read_excel(excel_path, header=None)
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        sys.exit(1)

    layout_map = {}
    
    # Rows A-H
    row_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    
    # If we have a grid that looks like data (rows <= 8, cols <= 12)
    # and doesn't have explicit headers "A", "1", etc.
    # We assume top-left is A01.
    
    n_rows, n_cols = df.shape
    
    if n_rows <= 8 and n_cols <= 12:
        # Direct mapping
        for r in range(min(n_rows, 8)):
            for c in range(min(n_cols, 12)):
                val = clean_name(df.iloc[r, c])
                if val.lower() == 'nan' or val == '':
                    continue # Empty well
                
                well_id = f"{row_labels[r]}{c+1:02d}"
                layout_map[well_id] = val
        return layout_map
        
    # If larger, fall back to header detection (omitted for brevity as previous logic was complex and failed)
    # But for this task, the file is known to be 6x12
    
    # Re-implement header detection if needed, but let's stick to the direct map if it fits bounds.
    # The previous complex logic failed because it expected headers.
    
    raise ValueError(f"Unexpected layout shape {df.shape}. Expected <= 8 rows and <= 12 columns.")

def get_well_from_filename(filename):
    # Common patterns: Specimen_001_A01_001.fcs, A01.fcs, etc.
    # Look for [A-H] followed by 1-2 digits, surrounded by non-alphanumeric or start/end
    # Regex: (?:^|[^A-Z0-9])([A-H])[-_]?(\d{1,2})(?:$|[^A-Z0-9])
    
    base = os.path.basename(filename)
    # matches A1, A01, A-01, A_01 embedded in string
    matches = re.findall(r"([A-H])[-_]?(\d{1,2})", base, re.IGNORECASE)
    
    # Filter matches that look like well IDs (cols 1-12)
    valid_wells = []
    for r, c in matches:
        c_int = int(c)
        if 1 <= c_int <= 12:
            valid_wells.append(f"{r.upper()}{c_int:02d}")
            
    if not valid_wells:
        return None
    
    # If multiple, usually the last one is the most specific (e.g. if filename has other codes)
    # But typically HTS filenames are consistent.
    # Let's take the one that matches the typical A01 pattern most closely.
    return valid_wells[-1] 

def find_channel_index(flow_data, channel_name):
    """Finds the index of the channel in the FlowData object."""
    # flow_data.channels is a dict: {index: {'PnN': 'Name', 'PnS': 'Desc', ...}}
    # Keys might be 'pnn', 'pns' (lowercase) or 'PnN', 'PnS' depending on version.
    
    for idx, info in flow_data.channels.items():
        # Try both casing styles
        pnn = info.get('PnN') or info.get('pnn') or ''
        pns = info.get('PnS') or info.get('pns') or ''
        
        if channel_name.lower() in pnn.lower() or channel_name.lower() in pns.lower():
            return idx
            
    return None

def get_text_keyword(flow_data, keys):
    for key in keys:
        if key in flow_data.text:
            return flow_data.text[key]
        for t in flow_data.text:
            if t.lower() == key.lower():
                return flow_data.text[t]
    return None


def find_channel_by_name(flow_data, names):
    for idx, meta in flow_data.channels.items():
        pnn = (meta.get("pnn") or meta.get("PnN") or "").lower()
        pns = (meta.get("pns") or meta.get("PnS") or "").lower()
        for name in names:
            if name.lower() in pnn or name.lower() in pns:
                return idx
    return None


def read_fcs_data(fcs_path, channel_name, doublet_threshold=None):
    """
    Reads an FCS file and returns the events for the specified channel.
    """
    try:
        flow_data = flowio.FlowData(fcs_path)
        idx = find_channel_index(flow_data, channel_name)
        
        if idx is None:
            # List available channels for debugging
            avail = []
            for i, meta in flow_data.channels.items():
                pnn = meta.get('PnN') or meta.get('pnn')
                pns = meta.get('PnS') or meta.get('pns')
                avail.append(f"{pnn}/{pns}")
            print(f"Warning: Channel '{channel_name}' not found in {os.path.basename(fcs_path)}. Available: {avail}")
            return None
            
        # events is a flat list, need to reshape
        events = np.array(flow_data.events) 
        # flowio events are usually flat, reshaped by channel count
        # But recent flowio might return structured data or flat.
        # check dimensions
        num_channels = len(flow_data.channels)
        num_events = len(events) // num_channels
        events = events.reshape((num_events, num_channels))
        
        # Retrieve data for channel (indices in flowio are 1-based in the dict keys usually, but array is 0-based?)
        # flowio channel keys are usually integer indices corresponding to the order in the file (1-based in standard, but check flowio docs)
        # In flowio 1.0+, keys are string indices '1', '2', etc. 
        # Let's verify the index.
        
        # Actually, flowio.FlowData.events is a 1D array.
        # We reshaped it. The channel index in the array corresponds to the order in the file.
        # The keys in `flow_data.channels` correspond to the `$Pn` keywords, usually '1', '2'...
        
        # If we found idx via the keys, we need to map key '1' to array index 0.
        array_idx = int(idx) - 1
        channel_data = events[:, array_idx]

        time_idx = find_channel_by_name(flow_data, ["Time"])
        time_data = None
        if time_idx:
            time_data = events[:, int(time_idx) - 1]

        timestep_value = get_text_keyword(flow_data, ["timestep", "$TIMESTEP"])
        try:
            timestep_value = float(timestep_value)
        except (TypeError, ValueError):
            timestep_value = None

        tot_value = get_text_keyword(flow_data, ["tot", "$TOT"])
        try:
            tot_value = int(float(tot_value))
        except (TypeError, ValueError):
            tot_value = len(channel_data)

        duration = None
        if time_data is not None and timestep_value:
            duration = (time_data.max() - time_data.min()) * timestep_value

        events_per_second = None
        if duration and duration > 0:
            events_per_second = tot_value / duration

        fsca_idx = find_channel_by_name(flow_data, ["FSC-A"])
        fsch_idx = find_channel_by_name(flow_data, ["FSC-H"])
        doublet_count = 0
        if fsca_idx and fsch_idx:
            fsca = events[:, int(fsca_idx) - 1]
            fsch = events[:, int(fsch_idx) - 1]
            ratio = fsca / (fsch + 1e-9)
            if doublet_threshold is not None:
                doublet_count = int(np.sum(ratio > doublet_threshold))

        return {
            "channel": channel_data,
            "time": time_data,
            "timestep": timestep_value,
            "tot": tot_value,
            "duration": duration,
            "events_per_second": events_per_second,
            "num_events": num_events,
            "doublet_count": doublet_count,
        }
        
    except Exception as e:
        print(f"Error reading {fcs_path}: {e}")
        return None

def biexponential_transform(events, a=0.5, b=1, c=0.5, d=1, f=0, w=0):
    """
    Applies Logicle (Biexponential) transform to events.
    We will use a simplified version: arcsinh(x / cofactor) which is equivalent to Logicle with specific params.
    Standard cofactor for Cytof/Flow is 150 or 5.
    Let's use cofactor=150 as a start.
    """
    # Simplified biexponential (arcsinh)
    # T(x) = asinh(x / cofactor)
    return np.arcsinh(events / 150.0)

def main():
    args = parse_arguments()
    args.wt_name   = clean_name(args.wt_name)
    args.blank_name = clean_name(args.blank_name)
    # Create output dir if not exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # 1. Load Plate Layout
    layout = read_plate_layout(args.layout_file)
    print(f"Loaded layout with {len(layout)} annotated wells.")
    print(f"Unique samples found: {set(layout.values())}")
    
    # 2. Find FCS files
    fcs_files = []
    if os.path.isdir(args.fcs_input):
        # scan directory for .fcs or .facs
        fcs_files = glob.glob(os.path.join(args.fcs_input, "*.fcs")) + glob.glob(os.path.join(args.fcs_input, "*.facs"))
    elif os.path.isfile(args.fcs_input):
        fcs_files = [args.fcs_input]
    else:
        print("Invalid FCS input path.")
        sys.exit(1)
        
    print(f"Found {len(fcs_files)} FCS files.")
    
    # 3. Process Data
    results = [] # list of dicts: {Well, Sample, Mean, SD, Events}
    
    for fpath in fcs_files:
        well = get_well_from_filename(fpath)
        if not well:
            continue
            
        if well not in layout:
            continue # Skip wells not in layout (or implied empty)
            
        sample = clean_name(layout[well])
        if str(sample).upper() == 'BLANK' and args.blank_name.upper() != 'BLANK':
             # If the layout says BLANK literally, but user provided a specific blank name, treat as blank? 
             # Or maybe user means "BLANK" in the layout IS the blank. 
             # The prompt says: "annotated with name or BLANK that means not measured".
             # But also "The fourth input is the name of the blank to use to correct for background activity".
             # This implies one of the "names" is the blank sample.
             # AND "BLANK" means "not measured" (empty well).
             # So we should SKIP wells named "BLANK".
             pass
        
        if str(sample).upper() == 'BLANK':
            continue # Not measured
            
        events_info = read_fcs_data(fpath, args.channel, args.doublet_threshold)
        if events_info is None:
            continue
        channel_events = events_info['channel']
        if len(channel_events) == 0:
            continue
        median_val = np.median(channel_events)
        sd_val = np.std(channel_events)
        total_events = events_info['tot']
        doublet_count = events_info.get('doublet_count', 0)
        corrected_events = total_events + doublet_count
        duration = events_info.get('duration')
        events_per_second = events_info.get('events_per_second')
        
        sampled_volume = None
        cells_in_well = None
        if args.flow_rate and duration and duration > 0:
            sampled_volume = args.flow_rate * duration
            if sampled_volume > 0:
                concentration = corrected_events / sampled_volume
                cells_in_well = concentration * args.well_volume

        results.append({
            'Well': well,
            'Sample': sample,
            'Mean': median_val, # Storing Median in the 'Mean' column to minimize refactoring, or rename?
                                # Let's rename to minimize confusion but update references.
            'Median': median_val,
            'SD_dist': sd_val, # SD of the distribution in the well
            'Events': channel_events,
            'Events_per_second': events_per_second,
            'Total_Events': total_events,
            'Corrected_Events': corrected_events,
            'Duration': duration,
            'Sampled_Volume_uL': sampled_volume,
            'Doublets': doublet_count,
            'Cells_in_Well': cells_in_well if cells_in_well is not None else np.nan
        })
        
    if not results:
        print("No matching data found.")
        sys.exit(1)
        
    df_res = pd.DataFrame(results)
    
    # 4. Identify Blank and WT
    # Filter for Blank samples
    blank_samples = df_res[df_res['Sample'] == args.blank_name]
    if blank_samples.empty:
        print(f"Error: No samples found matching blank name '{args.blank_name}'.")
        sys.exit(1)
        
    background_val = blank_samples['Median'].mean() # Average of the Medians of blanks
    print(f"Background Median ({args.blank_name}): {background_val:.2f}")
    
    # Filter for WT
    wt_samples = df_res[df_res['Sample'] == args.wt_name]
    if wt_samples.empty:
        print(f"Error: No samples found matching WT name '{args.wt_name}'.")
        sys.exit(1)
        
    # 5. Apply Correction
    df_res['Corrected_Val'] = df_res['Median'] - background_val
    
    # Calculate WT Value Corrected (Mean of WT replicates medians)
    wt_corrected_val = df_res[df_res['Sample'] == args.wt_name]['Corrected_Val'].mean()
    print(f"WT Corrected Median: {wt_corrected_val:.2f}")
    
    # 6. Group by Sample for Bar Charts (Triplicates)
    # We want grouping by Sample -> Mean of Means, SD of Means
    grouped = df_res.groupby('Sample').agg(
        Mean_of_Medians=('Median', 'mean'),
        SD_of_Medians=('Median', 'std'), 
        Corrected_Mean_of_Medians=('Corrected_Val', 'mean')
    ).reset_index()
    
    # Calculate Inverse FC
    # "Inverse fold change to the WT"
    # If WT = 100, Sample = 50. FC (Sample/WT) = 0.5. Inverse (WT/Sample) = 2.
    # We will use WT / Sample as requested.
    # Calculate for EACH replicate to get error bars.
    
    df_res['Inverse_FC'] = wt_corrected_val / df_res['Corrected_Val']
    
    # Add Replicate Number
    df_res['Rep_Num'] = df_res.groupby('Sample').cumcount() + 1

    # PLOTTING
    sns.set_theme(style="whitegrid")
    
    # Plot 1: Raw Measurements (Median +/- SD of Medians)
    plt.figure(figsize=(12, 6))
    
    # Sorting Logic for Raw Measurements:
    # 1. "Rep" and "WT" first (leftmost)
    # 2. Rest sorted Ascending by Median
    
    # Create a sorting key
    # We want Rep, then WT, then others sorted by value.
    # Let's assign explicit ranks: Rep=0, WT=1, Others=2.
    # Then sort by (Rank, Median).
    
    # Calculate mean median per sample for sorting
    sample_medians = df_res.groupby('Sample')['Median'].mean().reset_index()
    
    def get_rank(name):
        if name == args.blank_name: return 0
        if name == args.wt_name: return 1
        return 2
        
    sample_medians['Rank'] = sample_medians['Sample'].apply(get_rank)
    
    # Sort: First by Rank (Ascending), then by Median (Ascending)
    sample_medians = sample_medians.sort_values(['Rank', 'Median'], ascending=[True, True])
    raw_order = sample_medians['Sample'].tolist()
    
    # Define colors: Orange for Rep/WT, Blue (default) for others
    # We need a color mapping or list of colors matching the order
    # seaborn palette 'tab10': Blue is 0, Orange is 1.
    # We want Rep/WT -> Orange, Rest -> Blue? Or specific colors?
    # "Rep and WT in orange".
    
    colors = []
    for name in raw_order:
        if name in [args.blank_name, args.wt_name]:
            colors.append('orange') # Or 'tab:orange'
        else:
            colors.append('steelblue') # Or 'tab:blue'
            
    sns.barplot(data=df_res, x='Sample', y='Median', errorbar='sd', capsize=.1, order=raw_order, palette=colors, alpha=0.7)
    # Add strip plot for individual dots
    sns.stripplot(data=df_res, x='Sample', y='Median', order=raw_order, color='black', size=4, jitter=True, alpha=0.6)
    
    plt.title(f"Raw Measurements ({args.channel}) - Median +/- SD")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "1_raw_measurements.png"), dpi=300)
    plt.close()
    
    # Plot 2: Inverse Fold Change
    plt.figure(figsize=(12, 6))
    
    # Filter out Blank sample for Inverse FC plot to avoid infinity/scaling issues
    df_plot_ifc = df_res[df_res['Sample'] != args.blank_name].copy()
    
    # Debug print
    print("\nDEBUG: Individual Data for Inverse FC (Top 10):")
    print(df_plot_ifc[['Sample', 'Corrected_Val', 'Inverse_FC']].head(10))
    
    # Sort samples by the mean Inverse FC (descending order - "how good mutants are")
    # Assuming higher Inverse FC = "better" (more inhibition/lower signal than WT if WT is high?)
    # Or if standard FC, higher is better? 
    # "Inverse fold change to the WT" -> WT / Sample.
    # If Sample is low (high inhibition?), Inverse FC is High.
    # So sorting descending puts the strongest inhibitors first.
    
    # Calculate mean per sample for sorting
    sample_order = df_plot_ifc.groupby('Sample')['Inverse_FC'].mean().sort_values(ascending=False).index.tolist()
    
    # Plot using the individual replicate values, allowing seaborn to calculate Mean and SD
    sns.barplot(data=df_plot_ifc, x='Sample', y='Inverse_FC', errorbar='sd', capsize=.1, order=sample_order, alpha=0.7)
    # Add strip plot for individual dots
    sns.stripplot(data=df_plot_ifc, x='Sample', y='Inverse_FC', order=sample_order, color='black', size=4, jitter=True, alpha=0.6)
    
    plt.title(f"Inverse Fold Change to WT (Corrected by {args.blank_name}) - Median Based")
    plt.ylabel("Inverse FC (WT / Sample)")
    plt.xticks(rotation=45, ha='right')
    plt.axhline(1, color='r', linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "2_inverse_fold_change.png"), dpi=300)
    plt.close()

    # Plot 6: Inverse FC - Top 5 Mutants (bars + dots)
    top5 = sample_order[:5]
    plt.figure(figsize=(12, 6))
    group_df = df_plot_ifc[df_plot_ifc['Sample'].isin(top5)]
    order_top5 = [s for s in top5 if s in group_df['Sample'].unique()]
    sns.barplot(data=group_df, x='Sample', y='Inverse_FC', errorbar='sd', capsize=.1,
                order=order_top5, palette='tab20', alpha=0.8)
    sns.stripplot(data=group_df, x='Sample', y='Inverse_FC', order=order_top5,
                  color='black', size=4, jitter=True, alpha=0.7)
    plt.title("Top 5 Mutants by Inverse Fold Change")
    plt.axhline(1, color='r', linestyle='--')
    plt.ylabel("Inverse FC")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "6_inverse_fc_top5.png"), dpi=300)
    plt.close()

    # Plot 7: Inverse FC - Remaining Mutants
    rest = sample_order[5:]
    plt.figure(figsize=(14, 6))
    group_df = df_plot_ifc[df_plot_ifc['Sample'].isin(rest)]
    order_rest = [s for s in rest if s in group_df['Sample'].unique()]
    sns.barplot(data=group_df, x='Sample', y='Inverse_FC', errorbar='sd', capsize=.1,
                order=order_rest, palette='tab20', alpha=0.8)
    sns.stripplot(data=group_df, x='Sample', y='Inverse_FC', order=order_rest,
                  color='black', size=3, jitter=True, alpha=0.6)
    plt.title("Remaining Mutants (Inverse Fold Change)")
    plt.axhline(1, color='r', linestyle='--')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Inverse FC")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "7_inverse_fc_rest.png"), dpi=300)
    plt.close()

    # Additional plots: raw lowest 5 vs rest and fold change (sample/WT) lowest 5 vs rest
    raw_candidates = [s for s in raw_order if s not in [args.blank_name, args.wt_name]]
    raw_low5 = raw_candidates[:5]
    raw_rest = raw_candidates[5:]

    if wt_corrected_val == 0:
        df_res['Fold_Change'] = np.nan
    else:
        df_res['Fold_Change'] = df_res['Corrected_Val'] / wt_corrected_val

    def plot_group(group_samples, metric, title, filename):
        plt.figure(figsize=(12, 6))
        group_df = df_res[df_res['Sample'].isin(group_samples)]
        order = [name for name in group_samples if name in group_df['Sample'].unique()]
        sns.barplot(data=group_df, x='Sample', y=metric, errorbar='sd', capsize=.1,
                    order=order, palette='tab20', alpha=0.8)
        sns.stripplot(data=group_df, x='Sample', y=metric, order=order,
                      color='black', size=4, jitter=True, alpha=0.6)
        if metric == 'Fold_Change':
            plt.axhline(1, color='r', linestyle='--')
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(args.output, filename), dpi=300)
        plt.close()

    if raw_low5:
        plot_group(raw_low5, 'Median', "Lowest 5 Raw Median Samples", "8_raw_lowest5.png")
        plot_group(raw_low5, 'Fold_Change', "Lowest 5 Fold Change vs WT", "10_fc_lowest5.png")
    if raw_rest:
        plot_group(raw_rest, 'Median', "Remaining Raw Median Samples", "9_raw_rest.png")
        plot_group(raw_rest, 'Fold_Change', "Remaining Fold Change vs WT", "11_fc_rest.png")
    
    # Plot 3: Histograms
    # "Histogram of the Blue-CA cell counts for every mutant each as a subplot in a larger plot"
    # Title: replicate number and name
    # This implies we plot ALL wells.
    # We need a grid.
    
    # Use strictly 8x12 grid for the histograms to match plate layout
    rows = 8
    cols = 12
    
    # Share x and y axes as requested
    fig, axes = plt.subplots(rows, cols, figsize=(30, 20), sharex=True, sharey=True)
    
    # Map well to axes
    # axes is 8x12 array
    
    # Sort by Well to ensure correct placement if iterating, 
    # but better to address by index.
    
    # Create a map of Well -> Data
    well_data_map = df_res.set_index('Well').to_dict('index')
    
    row_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    
    for r in range(rows):
        for c in range(cols):
            ax = axes[r, c]
            well_id = f"{row_labels[r]}{c+1:02d}"
            
            if well_id in well_data_map:
                row_data = well_data_map[well_id]
                events = row_data['Events']
                sample_name = row_data['Sample']
                rep_num = row_data['Rep_Num']
                
                # Plot histogram
                # "bioexponentially transform" as requested
                # We will transform the data before plotting
                
                events_trans = biexponential_transform(events)
                
                # Plot using seaborn histplot
                # Do NOT use log_scale=True since we manually transformed it
                sns.histplot(events_trans, ax=ax, bins=50, element="step", fill=True)
                
                ax.set_title(f"{sample_name}\n({well_id})", fontsize=6)
                
                # Only label edges
                ax.set_xlabel("Biexp" if r == rows-1 else "")
                ax.set_ylabel("Cnt" if c == 0 else "")
                
                # Set limits in transformed space
                # 100,000 raw -> arcsinh(100000/150) ~= arcsinh(666) ~= 7.2
                # 0 raw -> 0
                limit_trans = np.arcsinh(100000/150.0)
                ax.set_xlim(left=0, right=limit_trans)
                
                # Configure Ticks to show original values (10^2, 10^3, 10^4, 10^5)
                # We need FixedLocator at the transformed positions
                major_ticks_raw = [100, 1000, 10000, 100000]
                major_ticks_trans = np.arcsinh(np.array(major_ticks_raw) / 150.0)
                
                ax.xaxis.set_major_locator(ticker.FixedLocator(major_ticks_trans))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(["$10^2$", "$10^3$", "$10^4$", "$10^5$"]))
                
                # Minor ticks
                minor_ticks_raw = []
                for exp in [2, 3, 4]:
                    minor_ticks_raw.extend(np.arange(2, 10) * (10**exp))
                minor_ticks_raw = [x for x in minor_ticks_raw if x < 100000]
                minor_ticks_trans = np.arcsinh(np.array(minor_ticks_raw) / 150.0)
                
                ax.xaxis.set_minor_locator(ticker.FixedLocator(minor_ticks_trans))
                
                # Enable ticks visibility for all plots
                ax.tick_params(axis='both', which='both', labelbottom=True, labelleft=True, labelsize=6)
                ax.grid(True, which='minor', alpha=0.2)
                ax.grid(True, which='major', alpha=0.5)
                
            else:
                # Empty well
                ax.axis('off')
                # ax.text(0.5, 0.5, "Empty", ha='center', va='center', fontsize=6, transform=ax.transAxes)
                
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "3_histograms.png"), dpi=300)
    plt.close()
    
    # Plot 4: Heatmap of SD (Wacko identification)
    # We need to map SD_dist back to 8x12 grid
    plate_sd = np.zeros((8, 12))
    plate_sd[:] = np.nan
    
    # Prepare IQR heatmap as requested (5th plot)
    plate_iqr = np.zeros((8, 12))
    plate_iqr[:] = np.nan
    
    row_map = {c: i for i, c in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])}
    
    for idx, row_data in df_res.iterrows():
        w = row_data['Well']
        r_char = w[0]
        c_idx = int(w[1:]) - 1
        if r_char in row_map and 0 <= c_idx < 12:
            plate_sd[row_map[r_char], c_idx] = row_data['SD_dist']
            
            # Calculate IQR from events
            ev = row_data['Events']
            q75, q25 = np.percentile(ev, [75 ,25])
            iqr = q75 - q25
            plate_iqr[row_map[r_char], c_idx] = iqr
            
    plt.figure(figsize=(12, 8))
    sns.heatmap(plate_sd, annot=True, fmt=".1f", 
                xticklabels=[str(i) for i in range(1, 13)],
                yticklabels=list(row_map.keys()), cmap="viridis")
    plt.title(f"Heatmap of Intra-Well Standard Deviation ({args.channel})")
    plt.savefig(os.path.join(args.output, "4_sd_heatmap.png"), dpi=300)
    plt.close()
    
    # Plot 5: Heatmap of IQR
    plt.figure(figsize=(12, 8))
    
    # Cap IQR at 6500
    plate_iqr_capped = np.copy(plate_iqr)
    
    # Mask values > 6500 for special handling (blacked out)
    # We'll allow the heatmap to plot up to 6500, and values above will be masked or set to a specific color
    # Seaborn heatmap 'mask' argument hides cells.
    # Or we can set vmax=6500 and set the 'over' color.
    
    # Let's use the 'over' color approach with a custom cmap
    import copy
    cmap = copy.copy(plt.get_cmap("RdBu_r")) # Blue-Red scale (reversed so Red is high usually? Or standard RdBu: Red high, Blue low? 
                                           # "blue redscale" -> usually means Blue=Low, Red=High?
                                           # sns 'coolwarm' or 'RdBu_r' (Red-Blue reversed -> Red high).
                                           # Let's use 'coolwarm' or 'RdBu_r'.
    cmap.set_over('black')
    cmap.set_bad('white') # For NaNs
    
    # Using masked array for heatmap? 
    # Actually, if we just clip the data for display, we lose the "blacked out" info unless we specifically handle it.
    # Best way: Use 'mask' for values > 6500 and plot a black background? 
    # Or simply use the 'vmax' and 'extend' features of matplotlib if passing ax.pcolormesh, but seaborn wraps it.
    # Seaborn doesn't support 'set_over' easily directly in heatmap without tweak.
    
    # Alternative: Mask the high values, plot the heatmap, then plot the high values as black rectangles?
    # Or simpler: set values > 6500 to NaN in data passed to heatmap, plot heatmap with 'bad' color transparent, 
    # and have a black background?
    # Let's try:
    # 1. Set > 6500 to NaN in a copy.
    # 2. Create a mask where <= 6500 (True means masked).
    
    # Actually, seaborn heatmap `mask` argument: If True, data is not shown.
    # So we mask the VALID data? No, we mask the INVALID (black) data to show black underneath?
    # No, seaborn plots on an Axes. We can set the Axes facecolor to black.
    # Then mask the values > 6500.
    
    # Let's try the standard cmap 'over' method first, it usually works if we pass `vmax` and `cbar_kws`.
    
    mask_over = plate_iqr > 6500
    
    # We want "blacked out with the values in there".
    # So we want to SEE the number, but the cell background is black.
    
    # 1. Plot the heatmap with vmax=6500.
    #    This will saturate >6500 to the max color (Red).
    #    We want it Black.
    
    # 2. Custom logic:
    #    Set cells > 6500 to NaN in the plotting array.
    #    Set background of axes to Black.
    #    But NaN usually maps to 'white' or transparent.
    #    If transparent, black background shows through.
    
    plate_iqr_plot = np.copy(plate_iqr)
    plate_iqr_plot[plate_iqr > 6500] = np.nan
    
    ax = sns.heatmap(plate_iqr_plot, annot=plate_iqr, fmt=".1f", # Annotate with ORIGINAL values
                xticklabels=[str(i) for i in range(1, 13)],
                yticklabels=list(row_map.keys()), 
                cmap="coolwarm", # Blue-Red
                vmax=6500,
                cbar_kws={'extend': 'max'}) # Show that it extends
    
    # Now, for the NaNs (which were > 6500), we want them Black.
    # The heatmap leaves them empty (white/transparent).
    # We can color them manually or set the facecolor.
    ax.set_facecolor('black') # This makes ALL NaNs black, including empty wells.
    # We don't want empty wells (true NaNs) to be black?
    # User said "blacked out with the values in there". 
    # Empty wells don't have values.
    # So we need to distinguish "Empty Wells" vs "High IQR Wells".
    
    # Revert approach:
    # Create a discrete colormap? No, too complex.
    # Let's iterate and color the high cells black manually after plotting?
    # Or: 
    # 1. Plot heatmap.
    # 2. Overlay black rectangles on >6500 cells.
    # 3. Re-add text annotation on top.
    
    # Let's stick to the Axes facecolor approach but mask ONLY the high values.
    # True empty wells (NaN in original) should remain white?
    # If we set facecolor to black, all NaNs are black.
    # We can plot the True Empty wells as White rectangles?
    
    # Better:
    # Use a custom colormap where the 'bad' (NaN) value is black?
    # But we have two types of 'bad': >6500 and Empty.
    
    # Let's go with: Mask > 6500. Plot.
    # Then overlay black patches for > 6500.
    
    # Reset figure to ensure clean state
    plt.clf()
    plt.figure(figsize=(12, 8))
    
    # Plot the main heatmap with data clamped to 6500 for color scaling
    # But we want the 'over' cells to be distinct.
    
    # Let's just use the `mask` feature for > 6500.
    # `mask` takes a boolean array. True = Masked (not plotted).
    mask_high = (plate_iqr > 6500)
    
    # Plot the heatmap (masking high values)
    ax = sns.heatmap(plate_iqr, annot=True, fmt=".1f", 
                     xticklabels=[str(i) for i in range(1, 13)],
                     yticklabels=list(row_map.keys()), 
                     cmap="coolwarm", 
                     vmax=6500,
                     mask=mask_high)
    
    # Now fill the masked (high) values with black
    # We can use a second heatmap that is ALL black, masked by ~mask_high?
    # Or just iterate. 8x12 is small.
    
    # Get the mesh object from heatmap?
    # Simpler: Iterate and add patches.
    from matplotlib.patches import Rectangle
    
    for r in range(8):
        for c in range(12):
            val = plate_iqr[r, c]
            if not np.isnan(val) and val > 6500:
                # Add black rectangle
                # Heatmap coords: x=c, y=r (origin top-left usually in seaborn heatmap?)
                # Seaborn heatmap: x is col index + 0.5, y is row index + 0.5 for centers.
                # Rectangles: (c, r) with width 1, height 1.
                ax.add_patch(Rectangle((c, r), 1, 1, fill=True, color='black', edgecolor='none', zorder=1))
                
                # Re-add text (annot=True does it, but it might be hidden by patch or we need white text)
                # We need to add text manually on top.
                ax.text(c + 0.5, r + 0.5, f"{val:.1f}", 
                        color='white', ha='center', va='center', weight='bold', zorder=2)
            elif np.isnan(val):
                 # Ensure empty wells are white/empty (default)
                 pass

    plt.title(f"Heatmap of Intra-Well IQR ({args.channel}) - Capped at 6500")
    plt.savefig(os.path.join(args.output, "5_iqr_heatmap.png"), dpi=300)
    plt.close()

    # Plot 12: Heatmap of computed cells per well (requires flow rate)
    if args.flow_rate:
        plate_cells = np.full((8, 12), np.nan)
        for idx, row_data in df_res.iterrows():
            w = row_data['Well']
            r_char = w[0]
            c_idx = int(w[1:]) - 1
            if r_char in row_map and 0 <= c_idx < 12:
                cells = row_data.get('Cells_in_Well')
                if cells is not None and not np.isnan(cells):
                    plate_cells[row_map[r_char], c_idx] = cells

        plt.figure(figsize=(12, 8))
        sns.heatmap(plate_cells, annot=True, fmt=".0f", 
                    xticklabels=[str(i) for i in range(1, 13)],
                    yticklabels=list(row_map.keys()), cmap="magma")
        plt.title(f"Heatmap of Estimated Cell Count per Well (⚠️ ROUGH ESTIMATE ±50%)")
        plt.savefig(os.path.join(args.output, "12_cells_heatmap.png"), dpi=300)
        plt.close()

        # Plot 13: Heatmap of estimated OD600 (pre-dilution)
        # cells_in_well is total cells in well_volume µL (diluted)
        # Original cell_volume µL was diluted into well_volume µL total
        # Dilution factor = well_volume / cell_volume
        # concentration in original (cells/mL) = (cells / well_volume) * dilution_factor * 1000
        # OD600 = concentration / od_calibration
        
        DILUTION_FACTOR = args.well_volume / args.cell_volume
        
        plate_od = np.full((8, 12), np.nan)
        for idx, row_data in df_res.iterrows():
            w = row_data['Well']
            r_char = w[0]
            c_idx = int(w[1:]) - 1
            if r_char in row_map and 0 <= c_idx < 12:
                cells = row_data.get('Cells_in_Well')
                if cells is not None and not np.isnan(cells):
                    # cells is total in well_volume µL well
                    # concentration in well (cells/µL) = cells / well_volume
                    # concentration in original (cells/µL) = (cells / well_volume) * DILUTION_FACTOR
                    # concentration in original (cells/mL) = concentration_original_per_uL * 1000
                    conc_original_per_mL = (cells / args.well_volume) * DILUTION_FACTOR * 1000.0
                    od600 = conc_original_per_mL / args.od_calibration
                    plate_od[row_map[r_char], c_idx] = od600

        plt.figure(figsize=(12, 8))
        sns.heatmap(plate_od, annot=True, fmt=".2f", 
                    xticklabels=[str(i) for i in range(1, 13)],
                    yticklabels=list(row_map.keys()), cmap="YlOrRd")
        plt.title(f"Heatmap of Estimated OD600 Pre-Dilution (⚠️ ROUGH ESTIMATE ±50%)")
        plt.savefig(os.path.join(args.output, "13_od600_heatmap.png"), dpi=300)
        plt.close()

    # EXCEL EXPORT
    print("Generating summary Excel file...")
    
    # Prepare Data for Excel (Pivot Tables)
    # Rows: Replicates (1, 2, 3...)
    # Columns: Mutations (Sample names)
    
    # 1. Inverse Fold Change Pivot
    # df_res has 'Sample', 'Rep_Num', 'Inverse_FC'
    # We want Mutation as Header -> Columns=Sample
    # Replicates in same column -> This is ambiguous.
    # "headers are the mutations and each of the replicates are in the same column"
    # -> Column A: Mutation 1. Rows 2,3,4: Rep 1, 2, 3 values.
    # This corresponds to pivot(index=Rep_Num, columns=Sample).
    
    # Filter out blanks for Inverse FC sheet? Usually yes, or include everything.
    # Let's include everything available in df_res that has an Inverse_FC value.
    
    try:
        pivot_ifc = df_res.pivot(index='Rep_Num', columns='Sample', values='Inverse_FC')
        pivot_raw = df_res.pivot(index='Rep_Num', columns='Sample', values='Median')
        
        # Sort columns alphabetically
        pivot_ifc = pivot_ifc.sort_index(axis=1)
        pivot_raw = pivot_raw.sort_index(axis=1)
        
        excel_path = os.path.join(args.output, "summary.xlsx")
        csv_ifc_path = os.path.join(args.output, "summary_inverse_fc.csv")
        csv_raw_path = os.path.join(args.output, "summary_raw_median.csv")
        
        # Save as standard Excel
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            pivot_ifc.to_excel(writer, sheet_name='Inverse Fold Change')
            pivot_raw.to_excel(writer, sheet_name='Raw Measurements Median')
            
        # Save as CSVs with European format (semicolon separator, comma decimal)
        # Excel "European" CSVs usually use ; as separator and , as decimal
        pivot_ifc.to_csv(csv_ifc_path, sep=';', decimal=',')
        pivot_raw.to_csv(csv_raw_path, sep=';', decimal=',')
            
        print(f"Summary Excel saved to {excel_path}")
        print(f"European CSV summaries saved to {csv_ifc_path} and {csv_raw_path}")
        
    except Exception as e:
        print(f"Error creating Excel file: {e}")

    print(f"Analysis complete. Outputs saved to {args.output}")
    
    # Print warning about cell count estimation accuracy
    if args.flow_rate:
        print("\n" + "="*70)
        print("⚠️  WARNING: Cell Count / OD600 Estimation Accuracy")
        print("="*70)
        print("""
The cell counts and OD600 values in plots 12-13 are ROUGH ESTIMATES only.

The BD LSR Fortessa does NOT have volumetric counting capability.
Our calculation uses: cells = (events / (flow_rate × time)) × well_volume

Known sources of error (estimated total: ±30-50%):
  • Flow rate not recorded in FCS file - user-provided value may be inaccurate
  • Flow rate varies during acquisition (±10-20%)
  • Doublet discrimination threshold is approximate
  • OD600 calibration factor varies by strain/conditions (±2-5×)

For accurate absolute counts, use:
  1. Counting beads (TruCount, CountBright) - Gold standard
  2. Volumetric cytometer (BD Accuri C6)
  3. Validate with hemocytometer or plate reader
""")
        print("="*70 + "\n")

if __name__ == "__main__":
    main()

