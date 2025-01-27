# Land use gradients drive spatial variation in Lassa fever host communities in the Eastern Province of Sierra Leone.

## Overview

This repository contains the code and data used for the analysis presented in "Land use gradients drive spatial variation in Lassa fever host communities in the Eastern Province of Sierra Leone." The study examines the impact of land use and species interactions on rodent community composition and Lassa fever spillover risk using hierarchical occupancy models.

## Repository structure

📂 **land-use-lassa-hosts**  
├── 📂 **data**                     # Raw and processed data  
│   ├── 📂 **geodata**              # Geographic data files*
|   |   |── 📂  **gadm**            # Sierra Leone shapefiles*
|   |   |── 📄  **iucn_clr.csv**    # IUCN colour table
|   |   |── 📂  **pop**             # Gridded population of the world population raster*
|   |   └── 📂  **wc2.1_tiles**     # Worldclim raster*
│   ├── 📂 **input**                # Original data files  
|   |   |── 📄  **rodent_data.csv** # Individual level rodent data
|   |   └── 📄  **trap_data.csv**   # Row per trap data
│   ├── 📂 **processed_data**       # Pre-processed data for analysis 
|   |   |── 📄  **descriptive_data.rds** # Data object matching rodent and trap data*
|   |   |── 📄  **detection_covariates.rds** # Trap and grid level detection covariate data*
|   |   └── 📄  **occurrence_covariates.rds** # Trap and grid level occurrence covariate data*
│   └── 📂 **output**               # Output data files and models  
|   |   |── 📄  **consolidated_rodent_trap.csv** # Cleaned matched rodent trap data*
|   |   └── 📂  **model**           # Models
|   |   |   |── **final_model.rds** # Final model**
|   |   |   └── **model_comparison.rds**  # Table produced from model comparisons
├── 🛠️ **R**                        # Analysis scripts  
│   ├── 📄 **00_setup.R**           # Load packages and project wide definitions
│   ├── 📄 **01_load_data.R**       # Load input data and process model covariates
│   ├── 📄 **02_descriptive_data.R**# Descriptive data analysis and visualisation
│   ├── 📄 **03_model_species_occ.R** # Model script
│   ├── 📄 **04_interpretation.R**  # Model checks, interpretation and visualisation code
│   ├── 📄 **05_co-occurrence.R**   # Co-occurrence analysis
│   └── 📄 **06_figure_1_script.R** # Produce Figure 1
├── 📊 **output**                   # Outputs generated from analysis  
│   ├── 📂 **figures/**         # Plots, graphs and supplementary figures
│   └── 📂 **tables/**          # Summary tables  
├── 📂 **report**                # Supporting documentation  
│   ├── 📄 **Supporting_Information_1.Rmd** # Study protocol
│   ├── 📄 **citations.bib**    # Manuscript references
│   └── 📄 **blinded_submission_manuscript.Rmd**  # Submitted blinded manuscript
├── 📜 **LICENSE**              # License information  
├── 📄 **README.md**            # Project description (this file)  
└── 🛠️ **session_info.txt**     # Environment details for reproducibility  

* Produced by running through the scripts
** Available from the [OpenScienceFramework](https://osf.io/jbm6y/?view_only=8f32ac4d8659464ca468914b8f89ae99)


