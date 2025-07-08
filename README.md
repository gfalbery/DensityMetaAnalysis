# DensityMetaAnalysis

## R code and datasets for “Density-dependent network structuring within and across wild animal systems”.

This repository provides all necessary materials to reproduce the meta-analysis and visualizations from the study by Albery et al. (2025), which investigates how local population density affects individual-level spatial and social network connectivity across diverse wildlife systems.

## Overview
This project performs a meta-analysis of individual-year-level data from 36 wildlife systems across 30 species. It examines the effects of local density on network centrality using both linear mixed models and generalized additive mixed models (GAMMs). The analysis distinguishes between spatial and social networks to understand how different types of networks react differently to density.

## Attribution
If using this code or approach, please cite:

Albery, G.F., et al. (2025). Density-dependent network structuring within and across wild animal systems. bioRxiv. doi:10.1101/2024.06.28.601262


## Code structure

R/
├── 00_DanFunctions.R                     # Helper functions used for the meta-analysis
├── 00_Phylopic Setup.R                   # Downloads and prepares animal silhouettes for figures

├── 01_Summarising System Replicates.R   # Aggregates system-level data and replicates
├── 01b_Calculating Proportion Area Overlap.R   # Calculates spatial overlap metrics for density validation

├── 02a_Linear Models.R                   # Fits linear models to system-level data
├── 02b_GAMMs.R                           # Fits Generalized Additive Mixed Models
├── 02c_Inflection Models.R               # Fits different models to different inflection points

├── 03a_Summarising Linear Models.R      # Extracts and summarises linear model results
├── 03b_Summarising GAMMs.R              # Extracts and summarises GAMM outputs
├── 03c_Summarising Inflection Models.R  # Extracts and summarises inflection point analyses

├── 04_Combining Summary Data.R          # Merges outputs from all systems for meta-analysis

├── 05a_Meta-Analysis.R                  # Runs phylogenetic meta-analysis on linear model slopes
├── 05b_Saturation Meta-Analysis.R       # Meta-analysis of GAMM smooths for saturation effects
├── 05c_Full Saturation Meta-Analysis.R  # Extended saturation analysis across contact types
├── 05d_Strength Meta-Analysis.R         # Meta-analysis of network strength

├── 06a_Figures.R                        # Generates manuscript figures
├── 06b_Density Schematic.R              # Creates schematic figure of density-network relationships

├── 07_Results.R                         # Formats key statistics and effect sizes for reporting

├── 08_Heteroskedasticity Validations.R # Validates model assumptions around variance structure
