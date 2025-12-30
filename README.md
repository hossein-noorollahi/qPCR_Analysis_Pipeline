# Comprehensive qPCR Analysis Pipeline in R

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
*(Note: Link this repository to Zenodo to generate your unique DOI badge)*

## Overview
This repository hosts a robust, automated R pipeline designed for the rigorous analysis of **Quantitative Real-Time PCR (qPCR)** data. Tailored for molecular biology research, this tool streamlines the workflow from raw $C_t$ values to high-resolution, publication-ready figures.

The script relies on the **Livak ($2^{-\Delta\Delta C_t}$)** method for normalization and incorporates **Welch's t-test** (unequal variances) to ensure statistical reliability.

## Key Features
*   **Flexible Data Input:** Accepts raw data in standard CSV format (Demo data included).
*   **Statistical Rigor:** 
    *   Calculates **Fold Change** & **log2 Fold Change**.
    *   Performs **Welch's t-test** for statistical significance ($p$-values).
    *   Computes **95% Confidence Intervals (CI)** for error estimation.
*   **Publication-Quality Visualization:**
    *   **Global Bar Plots:** High-contrast profiles with error bars.
    *   **Volcano Plots:** Dynamic thresholding for Up/Down-regulated genes.
    *   **Heatmaps:** Expression gradients with significance markers.
    *   **Individual Gene Plots:** Auto-generated separate plots for detailed inspection.

## Getting Started

### Prerequisites
You need **R** (version 4.0+) and the following packages:
```r
install.packages(c("tidyverse", "ggrepel", "RColorBrewer", "scales", "gridExtra"))