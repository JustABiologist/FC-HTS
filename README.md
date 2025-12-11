# FC-HTS-Fortessa

High-Throughput Screening (HTS) analysis pipeline for BD Fortessa flow cytometry data.

## Description

This pipeline processes 96-well plate flow cytometry data (FCS files), performs background correction, and generates various visualizations to analyze screening results.

Key features:
- **Automated Plate Layout Parsing**: Reads 96-well plate layouts (including 6x12 partial plates) from Excel.
- **Biexponential Transformation**: Applies arcsinh transformation to flow cytometry data for visualization.
- **Statistical Analysis**: Calculates Median intensity, Standard Deviation, and Inverse Fold Change relative to Wild Type (WT).
- **Visualizations**:
    - Raw Measurements Bar Chart (Median +/- SD, with individual points)
    - Inverse Fold Change Bar Chart (Sorted by magnitude, with individual points)
    - 96-well Grid Histograms (Biexponential axis)
    - Heatmaps for Intra-well SD and IQR (Quality Control)

## Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/FC-HTS-Fortessa.git
    cd FC-HTS-Fortessa
    ```

2.  Create the environment using Conda (Recommended):
    ```bash
    conda env create -f environment.yaml
    conda activate hts_pipeline
    ```

    Or using Python venv:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    ```

## Usage

Run the `analyze_hts.py` script with the required arguments:

```bash
python analyze_hts.py <layout_file> <fcs_directory> <wt_name> <blank_name> [--channel CHANNEL_NAME] [--output OUTPUT_DIR]
```

### Arguments:
- `layout_file`: Path to the Excel file (.xlsx) containing the plate layout.
- `fcs_directory`: Path to the folder containing the .fcs files.
- `wt_name`: The name used for the Wild Type sample in the layout (e.g., "WT").
- `blank_name`: The name used for the Blank/Negative control in the layout (e.g., "Rep").
- `--channel`: (Optional) The channel name to analyze (default: "Blue C-A").
- `--flow-rate`: Flow rate in µL/sec used for converting events to concentration (required for cell/OD outputs).
- `--doublet-threshold`: FSC-A/FSC-H ratio above which an event is treated as a doublet and counted twice (default: 1.5).
- `--od-calibration`: Cells/mL per OD600 unit for your organism (default: 8e8 for E. coli in LB).
- `--cell-volume`: Volume of cell suspension added to each well in µL (default: 20).
- `--well-volume`: Total volume in each well in µL (default: 300).
- `--output` / `-o`: (Optional) Output directory for results (default: current directory).

### Example:

```bash
python analyze_hts.py "data/layout.xlsx" "data/experiment_001/" WT Rep --channel "Blue C-A" --output results/
```

## Outputs

The script generates the following in the output directory:
1.  **Plots**:
    - `1_raw_measurements.png`: Median fluorescence intensity.
    - `2_inverse_fold_change.png`: Inverse fold change (WT / Sample).
    - `3_histograms.png`: 8x12 grid of histograms.
    - `4_sd_heatmap.png`: Heatmap of standard deviation.
    - `5_iqr_heatmap.png`: Heatmap of IQR (outliers > 6500 blacked out).
    - `6–11_*.png`: Supplemental split plots for Top-5 vs rest (inverse FC, raw medians, fold change).
    - `12_cells_heatmap.png`: Estimated total cells per well (⚠️ rough estimate ±50%).
    - `13_od600_heatmap.png`: Estimated OD600 of the undiluted cell medium (⚠️ rough estimate ±50%).

> **⚠️ Warning**: Cell counts and OD600 estimates (plots 12-13) are approximations only. The BD LSR Fortessa does not have volumetric counting capability. The calculation relies on user-provided flow rate which is not recorded in the FCS file. For accurate absolute counts, use counting beads (TruCount/CountBright) or a volumetric cytometer.
2.  **Data**:
    - `summary.xlsx`: Excel file with pivoted data (Inverse FC and Raw Medians).
    - `summary_inverse_fc.csv`: European-formatted CSV (semicolon separator).
    - `summary_raw_median.csv`: European-formatted CSV (semicolon separator).
