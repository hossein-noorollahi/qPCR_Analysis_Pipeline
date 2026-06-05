---
title: 'Advanced qPCR Analysis Pipeline: A Reproducible R Framework for Automated Gene Expression Quantification and Visualization'
tags:
  - R
  - qPCR
  - gene expression
  - Livak method
  - bioinformatics
  - data visualization
  - reproducibility
authors:
  - name: Hossein Noorollahi
    orcid: 0000-0001-5685-8821
    affiliation: 1
    corresponding: true
  - name: Mitra Heydari Nasrabadi
    affiliation: 1
    corresponding: true
  - name: Zarrin Minuchehr
    affiliation: 2
  - name: Somaye Ehtesham
    affiliation: 1
  - name: Bijan Bambai
    affiliation: 2
affiliations:
  - index: 1
    name: Department of Biology, Parand Branch, Islamic Azad University, Parand, Iran
  - index: 2
    name: Department of Systems Biotechnology, National Institute of Genetic Engineering and Biotechnology, Tehran, Iran
date: 5 June 2026
bibliography: paper.bib
---

# Summary

Quantitative real-time PCR (RT-qPCR) is a cornerstone technique for quantifying gene expression, and the comparative $C_t$ (Livak, $2^{-\Delta\Delta C_t}$) method is the most widely used mathematical model for relative quantification [@Livak2001; @Bustin2009]. Despite the mathematical simplicity of this method, its practical implementation often relies on manual calculations in spreadsheet software, which is time-consuming, error-prone, and hampers computational reproducibility [@Ziemann2016].

The **Advanced qPCR Analysis Pipeline** is an open-source, fully automated R script that transforms raw $C_t$ data (standard CSV format with columns `Sample`, `Target`, `C_t`) into statistically validated results and publication-ready figures. The pipeline automatically identifies the reference gene and treatment groups, computes $\Delta C_t$ and $\Delta\Delta C_t$, and calculates relative fold changes using the Livak method. It then performs a robust statistical analysis using **Welch’s two-sample t-test** (which does not assume equal variances) to derive exact p-values, standard errors, and 95% confidence intervals. Finally, it generates high-resolution visualizations using the `ggplot2` and `ggrepel` R packages: global bar plots with error bars, volcano plots with collision-free labels, and heatmaps annotated with significance asterisks.

# Statement of need

Manual processing of qPCR data in spreadsheets remains common practice, yet it introduces several critical problems. First, copy-and-paste errors and incorrect formula propagation are frequent and often remain undetected [@Ziemann2016]. Second, researchers routinely apply Student’s t-test without checking the assumption of equal variances, which can inflate the type I error rate when biological variances differ between control and treatment groups. Third, generating publication-quality graphics typically requires moving between multiple software tools (e.g., spreadsheet, graphing software, illustration program), increasing the risk of inconsistencies and slowing down research.

The **Advanced qPCR Analysis Pipeline** directly addresses these issues by providing a single, script-based, reproducible solution. It enforces statistical rigor by using **Welch’s t-test** [@Ruxton2006], which is more appropriate when variances are unequal. It automatically computes standard errors and confidence intervals, and outputs a comprehensive CSV report containing fold change, $\log_2$ fold change, p-values, and regulatory status for each gene. Moreover, it seamlessly integrates statistical testing with visualization, producing high-resolution plots ready for publication without any manual editing.

Existing R packages for qPCR analysis often require manual coding for each new dataset or produce only limited visual outputs. Our pipeline is designed to work out-of-the-box with a simple CSV file, making it accessible to molecular biologists with minimal R experience. The code is thoroughly commented and accompanied by a dummy dataset, allowing users to test the complete workflow instantly. By eliminating manual spreadsheet manipulation and integrating state-of-the-art statistical testing, this pipeline reduces analytical bottlenecks and enhances the reproducibility of gene expression studies. Crucially, by dynamically adjusting for inter-patient variance, this pipeline actively mitigates the risk of reporting false-positive clinical biomarkers in heterogeneous datasets.

# Authorship Contributions and Intellectual Property

**Hossein Noorollahi** independently conceived the original idea for this software, conceptualized its mathematical framework, and served as the sole architect and lead developer responsible for the entirety of the R programming. The co-authors acted in a supervisory and advisory capacity, providing academic consultation and oversight regarding the biological context of gene expression, without involvement in the software's ideation or coding. The software is distributed publicly under the open-source **MIT License**, with the intellectual and commercial rights of the source code retained solely by the lead developer (H.N.).

# Funding and Institutional Disclaimer

This software development project was conducted as an independent computational endeavor. The authors explicitly declare that no financial support, grants, or institutional resources were received for the ideation, development, or execution of this pipeline from any affiliated institutions, including the Islamic Azad University (Parand Branch) and the National Institute of Genetic Engineering and Biotechnology (NIGEB).

# Conflict of Interest

The authors declare that the research and software development were conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

# Acknowledgements

We thank the open-source R community for developing the `tidyverse`, `ggplot2`, and `ggrepel` packages which served as foundational dependencies for the graphical components of this framework.

# References
