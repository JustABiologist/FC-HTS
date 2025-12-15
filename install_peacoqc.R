#!/usr/bin/env Rscript
# PeacoQC Installation Script
# Run this after creating the conda environment to install PeacoQC and dependencies

message(paste(rep("=", 70), collapse = ""))
message("PeacoQC Installation")
message(paste(rep("=", 70), collapse = ""))

# Check if BiocManager is available (should be installed via conda)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install Bioconductor dependencies (flowWorkspace is required by PeacoQC)
bioc_deps <- c("flowCore", "flowWorkspace", "ComplexHeatmap")
message("\nInstalling Bioconductor dependencies: ", paste(bioc_deps, collapse = ", "))
message("(This may take several minutes...)\n")

for (pkg in bioc_deps) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("  Installing", pkg, "..."))
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
        message(paste("  ✓", pkg, "already installed"))
    }
}

# Install PeacoQC from GitHub
message("\nInstalling PeacoQC from GitHub (saeyslab/PeacoQC)...")

if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools...")
    install.packages("devtools", repos = "https://cloud.r-project.org")
}

devtools::install_github("saeyslab/PeacoQC", upgrade = "never", quiet = FALSE)

# Verify installation
message(paste("\n", paste(rep("=", 70), collapse = ""), sep = ""))
if (requireNamespace("PeacoQC", quietly = TRUE)) {
    message("PeacoQC installed successfully!")
    message(paste("  Version:", packageVersion("PeacoQC")))
    message("\nInstalled packages:")
    for (pkg in c("flowCore", "flowWorkspace", "ComplexHeatmap", "PeacoQC")) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            message(paste("  OK:", pkg, "-", as.character(packageVersion(pkg))))
        }
    }
} else {
    stop("PeacoQC installation failed!")
}
message(paste(rep("=", 70), collapse = ""))
