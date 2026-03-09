# Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis

## Abstract 

Characterizing how biodiversity has changed through Earth’s history and uncovering the processes that have driven those changes remain a significant challenge. Haar fluctuation analysis, a recently developed time series method, has been suggested as a powerful new tool to infer macroevolutionary drivers and assess system stability. Yet the ability of this method to identify unique drivers or the timing and dynamics of biodiversity, particularly in biased time series, has not been demonstrated. Here, we assess Haar fluctuation analysis and cross-Haar correlations using process-based ecological simulations that incorporate realistic sampling and depositional biases. We find that simpler (neutral) mechanisms can produce patterns observed in the Phanerozoic record, and that uneven sampling and sedimentary hiatuses can distort scaling relationships, cautioning against mechanistic interpretations. Nonetheless, Haar fluctuation analysis can reliably distinguish stabilizing from non-stabilizing dynamics, even under severe sampling bias, supporting the identification of a long-term equilibrium in Phanerozoic marine biodiversity. Our results suggest that this aspect of Haar fluctuation analysis will be robust whenever the system’s return time exceeds the temporal grain of sampling, underscoring the broader value of timescale-based approaches for studying biodiversity dynamics. 

## Publication Information 

Petryshen, W., Pincelli, HM., Vasseur, D. 2026. Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis. Proc. R. Soc. B. DOI: 10.1098/rspb.2025.2840

## General Information

1. Datasets:
Each of the listed datasets is used either within the main article or supplementary information. To recreate the figures and analysis run the analysis or supplement Jupyter Notebooks.
    - Fig2_Data
    - Fig4_RW_Data
    - Fig4_UNTB_Data
    - Neotoma_Results
    - Plantonic_Results

2. Author Information:
	- Principal Investigator Contact Information
		Name: Wyatt Petryshen;
		Institution: Department of Earth and Planetary Sciences, Yale University;
		Email: wyatt.petryshen@yale.edu;

3. License information: 
    - CC-BY 4.0

4.  DOI for data and code repository:
    - Code Repository: Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis. 2025. Zenodo. https://doi.org/10.5281/zenodo.18759980

## Data and File Overview

This repository contains the code needed to perform the simulation experiments and analysis described in the above paper.

Please make sure all directory locations and save folders are correctly set within each notebook prior to running analysis. The repository also contains a requirements.txt file to ensure all analysis dependencies can be recreated. Data needs to be unzipped prior to running analysis.

### Run Simulations 
Simulations can be re-run using the files listed within the cluster folder. Please make sure to correctly set path variables and save folders. Simulation outputs should be saved to the same folder structure as contained within the data folder. These simulation results will then be automatically imported when running the Jupyter Notebooks within the analysis folder. 

### Data
If you choose not to rerun each simulation as described above, you can import all of the required data from within the data folder. 
- Data for Figure 2: 
    - Fig2_Data/RW: Required .npz files for RW model in figure 2 (imported in Analysis/Fig2_Analysis.ipynb)
    - Fig2_Data/UNTB: Required .npz files for UNTB model in figure 2 (imported in Analysis/Fig2_Analysis.ipynb)
- Data for Figure 4:
    - Fig4_UNTB_Data
        - Fig4_UNTB_Data/UNTB_Dep_Cut: Required .npz files for UNTB models in figure 4 (imported in Analysis/Fig4_Analysis.ipynb)
        - Ideal_Case_UNTB_Dep_Cut.npz: Required .npz files for the ideal UNTB case in figure 4 (imported in Analysis/Fig4_Analysis.ipynb)
    - Fig4_UNTB_Data
        - Fig4_UNTB_Data/RW_Dep_Cut: Required .npz files for RW models in figure 4 (imported in Analysis/Fig4_Analysis.ipynb)
        - Ideal_Case_RW_Dep_Cut.npz: Required .npz files for the ideal RW case in figure 4 (imported in Analysis/Fig4_Analysis.ipynb)
- Data for supplementary empirical analysis (SFig 14 - 19): 
    - Neotoma_Results: Required .json files for empirical analysis of tree pollen assemblages (imported in supplements/BioDeepTime_Analysis.ipynb)
    - Plantonic_Results: Required .json files for empirical analysis of planktonic foraminifera assemblages (imported in supplements/BioDeepTime_Analysis.ipynb)

Note: Column names for dataframes are included on data import for both the .json and .npz files. 

### Code 

#### Main Figures and Analysis
1. Analysis for Figure 2:
    - analysis/Fig2_Analysis.ipynb
2. Analysis for Figure 3:
    - analysis/Fig3_Analysis.ipynb
3. Analysis for Figure 4:
    - analysis/Fig4_Analysis.ipynb (depends on code within cluster folder)

#### Supplementary Figures and Analysis
0. Source code and user defined functions:
    - src 
1. Analysis for SFig 1 - 4: 
    - supplements/UNTB_Parameter_Choice.ipynb: Interactively create supplementary figures 1 - 4 using this Juypter Notebook; Parameters will need to be changed dependent upon which parameters are to be explored.
2. Analysis for SFig 5: 
    - analysis/Fig3_Analysis.ipynb
3. Analysis for SFig 6 - 8: 
    - analysis/Fig4_Analysis.ipynb
4. Analysis for SFig 9: 
    - supplements/Exogenous_Factor_Example.ipynb
5. Analysis for SFig 10 - 12: 
    - supplements/Depositional_Model.ipynb
6. Analysis for SFig 13: 
    - supplements/Haar_Algorithm_Comparison.ipynb
7. Analysis for SFig 14 - 19: 
    - supplements/BioDeepTime_Analysis.ipynb: Plotting of analysis
    - supplements/BioDeepTime_Neotoma.R: Fetching, cleaning, saving data
    - supplements/BioDeepTime_Planktonic_Forams.R: Fetching, cleaning, saving data
    - supplements/BioDeepTime_R_Analysis.R: Haar fluctuation analysis of data

## Methodological Information

This repository is organized such that each publication figure can be recreated within individual Juypter notebooks. All analysis can also be reconducted by running files in the cluster folder. All of the functions source code is contained within the src directory. The analysis conducted for Figure 4 was completed on an HPC system. As such, we have included the simulation output files used to construct the final figure within the data folder. Original simulation files for the HPC are included within the "cluster" folder.

To recreate the figures or analysis in the main article or supplements, simply run each Jupyter Notebook and its dependencies as described above. 

### Instrument- or software-specific information (including package versions) needed to interpret the data: 

### R
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20

- chronosphere 0.6.1
- ggrepel 0.9.6
- tidyverse 2.0.0
- dplyr 1.1.4
- zoo 1.8-14
- parallel 4.5.1
- stringr 1.5.2
- jsonlite 2.0.0
- divDyn 0.8.3
- ggplot2 4.0.1
- maps 3.4.3
- viridis 0.6.5
- ncdf4 1.24
- scales 1.4.0
- reshape2 1.4.4

### Python 
Python version: 3.13.9

All required python packages are contained within the requirements.txt file

