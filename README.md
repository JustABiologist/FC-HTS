# FC-HTS-Fortessa

High-Throughput Screening (HTS) analysis pipeline for BD Fortessa flow cytometry data.

## Description

This pipeline processes 96-well plate flow cytometry data (FCS files), performs background correction, and generates various visualizations to analyze screening results.

Key features:
- **Automated Plate Layout Parsing**: Reads 96-well plate layouts (including 6x12 partial plates) from Excel.
- **Biexponential Transformation**: Applies arcsinh transformation to flow cytometry data for visualization.
- **Statistical Analysis**: Calculates Median intensity, Standard Deviation, and Inverse Fold Change relative to Wild Type (WT).
- **Visualizations**:
    - Raw Measurements Bar Chart (Median +/- SD)
    - Inverse Fold Change Bar Chart (Sorted by magnitude)
    - 96-well Grid Histograms (Log-scale/Biexponential axis)
    - Heatmaps for Intra-well SD and IQR (Quality Control)

## Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/FC-HTS-Fortessa.git
    cd FC-HTS-Fortessa
    ```

2.  Create a virtual environment (recommended):
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Run the `analyze_hts.py` script with the required arguments:

```bash
python analyze_hts.py <layout_file> <fcs_directory> <wt_name> <blank_name> [--channel CHANNEL_NAME]
```

### Arguments:
- `layout_file`: Path to the Excel file (.xlsx) containing the plate layout.
- `fcs_directory`: Path to the folder containing the .fcs files.
- `wt_name`: The name used for the Wild Type sample in the layout (e.g., "WT").
- `blank_name`: The name used for the Blank/Negative control in the layout (e.g., "Rep").
- `--channel`: (Optional) The channel name to analyze (default: "Blue-CA").

### Example:

```bash
python analyze_hts.py "data/layout.xlsx" "data/experiment_001/" WT Rep --channel "Blue C-A"
```

## Outputs

The script generates the following plots in the current directory:
1.  `1_raw_measurements.png`: Bar chart of median fluorescence intensity.
2.  `2_inverse_fold_change.png`: Bar chart of inverse fold change (WT / Sample).
3.  `3_histograms.png`: 8x12 grid of histograms for each well.
4.  `4_sd_heatmap.png`: Heatmap of standard deviation within each well.
5.  `5_iqr_heatmap.png`: Heatmap of interquartile range (IQR) within each well (QC tool).

