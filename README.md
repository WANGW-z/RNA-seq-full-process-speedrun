# RNA-Seq Data Analysis Pipeline

A robust, scalable, and reproducible Snakemake pipeline for processing bulk RNA-Sequencing data from raw FASTQ files to differential expression and functional enrichment analysis.

## Overview

This repository contains an automated bioinformatics pipeline for comprehensive RNA-Seq analysis. It performs quality control, adapter trimming, alignment, gene quantification, differential expression testing, and extensive functional enrichment (GSEA, GO, KEGG) with publication-ready visualizations. The pipeline is built using [Snakemake](https://snakemake.github.io/), ensuring reproducibility and ease of execution on high-performance computing clusters and workstations.

## Features

- **Full Automation**: Handles the entire analysis from raw data to results with a single command.
- **Reproducibility**: Managed Conda environments ensure consistent software versions.
- **Comprehensive QC**: Integrates FastQC and MultiQC for thorough quality assessment before and after trimming.
- **Standardized Analysis**: Utilizes established tools (HISAT2, featureCounts, DESeq2, clusterProfiler) for reliable results.
- **Rich Visualization**: Automatically generates PCA plots, heatmaps, volcano plots, and enrichment bar charts.
- **Flexible Configuration**: Easily adaptable to different projects via a central JSON config file.

## Workflow Overview

The pipeline performs the following steps:
1. **Quality Control (Raw Data)**: `FastQC` on raw FASTQ files.
2. **Trimming**: `Trimmomatic` to remove adapters and low-quality bases.
3. **Quality Control (Trimmed Data)**: `FastQC` on the cleaned FASTQ files.
4. **QC Report Aggregation**: `MultiQC` summarizes all QC results into a single report.
5. **Alignment**: `HISAT2` aligns reads to a reference genome.
6. **Quantification**: `featureCounts` counts reads mapping to genes.
7. **Differential Expression**: `DESeq2` identifies statistically significant differentially expressed genes.
8. **Functional Enrichment**: `clusterProfiler` performs GSEA, GO, and KEGG analysis on the results.

## Citation and Acknowledgement
If you use this pipeline in your research, please kindly include an acknowledgment in your publication:
We gratefully acknowledge the use of the RNA-Seq analysis pipeline provided by WEN WANG, Ph.D. (Beijing Institute of Ophthalmology, Beijing Tongren Hospital).
Researcher Information:
WANG, Wen, Ph.D.
Bioinformatic Researcher
Beijing Institute of Ophthalmology (BIO), Beijing Tongren Hospital, Capital Medical University
E-mail: 18519720962@163.com

## License
https://img.shields.io/badge/License-MIT-yellow.svg

This project is licensed under the MIT License - see the LICENSE.md file for details.

## Support
For questions, bugs, or suggestions regarding this pipeline, please open an [Issue](<your-github-repo-url>/issues) on GitHub or contact Dr. WEN WANG directly via email.

## Disclaimer: 
This pipeline is provided for research purposes only. Users are responsible for validating the results and ensuring the chosen parameters and reference files are appropriate for their specific experiment.
