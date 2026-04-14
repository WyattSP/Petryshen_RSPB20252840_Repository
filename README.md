# Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis

## Abstract 

Characterizing how biodiversity has changed through Earth’s history and uncovering the processes that have driven those changes remain a significant challenge. Haar fluctuation analysis, a recently developed time-series method, has been suggested as a powerful new tool to infer macroevolutionary drivers and assess system stability. Yet the ability of this method to identify unique drivers or the timing and dynamics of biodiversity, particularly in biased time series, has not been demonstrated. Here, we assess Haar fluctuation analysis and cross-Haar correlations using process-based ecological simulations that incorporate realistic sampling and depositional biases. We find that simpler (neutral) mechanisms can produce patterns observed in the Phanerozoic record, and that uneven sampling and sedimentary hiatuses can distort scaling relationships, cautioning against mechanistic interpretations. Nonetheless, Haar fluctuation analysis can reliably distinguish stabilizing from non-stabilizing dynamics, even under severe sampling bias, supporting the identification of a long-term equilibrium in Phanerozoic marine biodiversity. Our results suggest that Haar fluctuation analysis will be robust for detecting stability whenever the time series is of sufficient resolution relative to duration, and duration relative to the system’s return time. More broadly, these findings underscore the value of time-scale-based approaches for studying biodiversity dynamics.


## Publication Information 

Petryshen W, Hull PM, Vasseur DA. 2026 Detecting stabilizing dynamics in biased biodiversity time series using Haar fluctuation analysis. Proc. R. Soc. B 293: 20252840. https://doi.org/10.1098/rspb.2025.2840

## General Information

1. Datasets:
Each of the listed datasets are used either within the main article or supplementary information. To recreate the figures and analysis run either the analysis or supplemental Jupyter Notebooks. Although each dataset is provided, they can be created by running files in the cluster folder, or the BioDeepTime R scripts. Note that "allow_pickle=True" when importing .npz files with np.load(file_path, allow_pickle=True).
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

5. General Analysis Overview (Input -> Code -> Output):
    - Fig2_Data -> Fig2_Analysis.ipynb -> Figure 2
    - Fig3_Analysis.ipynb -> Figure 3; SFig 5
    - Fig4_RW_Data, Fig4_UNTB_Data -> Fig4_Analysis.ipynb -> Figure 4, SFig 6 - 8
    - UNTB_Parameter_Choice.ipynb -> SFig 1 - 4
    - Exogenous_Factor_Example.ipynb -> SFig 9
    - Depositional_Model.ipynb -> SFig 10 - 12
    - Haar_Algorithm_Comparison.ipynb -> SFig 13
    - Neotome_Results, Plantonic_Results -> BioDeepTime_Analysis.ipynb -> SFig 14 - 19

## Data and File Overview

This repository contains the code needed to perform the simulation experiments and analysis described in the above paper.

Please make sure all directory locations and save folders are correctly set within each notebook prior to running analysis. The repository also contains a requirements.txt file to ensure all analysis dependencies can be recreated. Data needs to be unzipped prior to running analysis.

### Run Simulations 
Simulations can be re-run using the files listed within the cluster folder. Please make sure to correctly set path variables and save folders. Simulation outputs should be saved to the same folder structure as contained within the data folder. These simulation results will then be automatically imported when running the Jupyter Notebooks within the analysis folder. 

### Data
If you choose not to rerun each simulation as described above, you can import all of the required data from within the data folder. 
- Data for Figure 2: 
    - Fig2_Data/RW: Required .npz files for RW model in figure 2 (imported by Analysis/Fig2_Analysis.ipynb)
    - Fig2_Data/UNTB: Required .npz files for UNTB model in figure 2 (imported by Analysis/Fig2_Analysis.ipynb)
- Data for Figure 4:
    - Fig4_UNTB_Data
        - Fig4_UNTB_Data/UNTB_Dep_Cut: Required .npz files for UNTB models in figure 4 (imported by Analysis/Fig4_Analysis.ipynb)
        - Ideal_Case_UNTB_Dep_Cut.npz: Required .npz files for the ideal UNTB case in figure 4 (imported by Analysis/Fig4_Analysis.ipynb)
    - Fig4_UNTB_Data
        - Fig4_UNTB_Data/RW_Dep_Cut: Required .npz files for RW models in figure 4 (imported by Analysis/Fig4_Analysis.ipynb)
        - Ideal_Case_RW_Dep_Cut.npz: Required .npz files for the ideal RW case in figure 4 (imported by Analysis/Fig4_Analysis.ipynb)
- Data for supplementary empirical analysis (SFig 14 - 19): 
    - Neotoma_Results: Required .json files for empirical analysis of tree pollen assemblages (imported by supplements/BioDeepTime_Analysis.ipynb)
    - Plantonic_Results: Required .json files for empirical analysis of planktonic foraminifera assemblages (imported by supplements/BioDeepTime_Analysis.ipynb)

Note: Column names for dataframes are included on data import for both the .json and .npz files. 

### Code 

#### Main Figures and Analysis
1. Analysis for Figure 2: Haar fluctuation analysis and cross-haar correlations for unbiased neutral and species-level models; all data is produced wihtin notebook or can be imported from the data folder
    - analysis/Fig2_Analysis.ipynb
2. Analysis for Figure 3: Haar fluctuation analysis and cross-haar correlations for neutral community after abundance bias; all data is produced within notebook
    - analysis/Fig3_Analysis.ipynb
3. Analysis for Figure 4: Haar fluctuation analysis and cross-haar correlations for neutral and species-level models under depositional bias; data can either be produced from files in the cluster folder, or imported from data folder
    - analysis/Fig4_Analysis.ipynb (depends on code within cluster folder)

#### Supplementary Figures and Analysis
0. Source code and user defined functions:
    - src 
1. Analysis for SFig 1 - 4: Explore parameter combinations for neutral and species-level models and bias; all required data is produced within notebook
    - supplements/UNTB_Parameter_Choice.ipynb: Interactively create supplementary figures 1 - 4 using this Juypter Notebook; Parameters will need to be changed dependent upon which parameters are to be explored.
2. Analysis for SFig 5: Moment convergence of community turnover distributions; contained at the end of the Fig3_Analysis.ipynb; all data is produced wihtin notebook
    - analysis/Fig3_Analysis.ipynb
3. Analysis for SFig 6 - 8: Exploration of community turnover distributions; contained at the end of the Fig4_Analysis.ipynb; data can either be produced from files in the cluster folder, or imported from data folder
    - analysis/Fig4_Analysis.ipynb
4. Analysis for SFig 9: Haar fluctuation analysis on a climatically driven community; all data is produced within notebook
    - supplements/Exogenous_Factor_Example.ipynb
5. Analysis for SFig 10 - 12: Comparisons of Pareto distributions, and Sadler plots for the depositional models; all required data is produced within notebook
    - supplements/Depositional_Model.ipynb
6. Analysis for SFig 13: Comparions of different Haar fluctuation analysis algorithms; all data is produced within notebook
    - supplements/Haar_Algorithm_Comparison.ipynb
7. Analysis for SFig 14 - 19: Emperical analysis of BioDeepTime data; the original dataset is imported from R via the chronosphere package; processed data also contained in data folder
    - supplements/BioDeepTime_Analysis.ipynb: Plotting of analysis
    - supplements/BioDeepTime_Neotoma.R: Fetching, cleaning, saving data
    - supplements/BioDeepTime_Planktonic_Forams.R: Fetching, cleaning, saving data
    - supplements/BioDeepTime_R_Analysis.R: Haar fluctuation analysis of data

## Methodological Information

This repository is organized such that each publication figure can be recreated within individual Juypter notebooks. All analysis is also reproducable by running files in the cluster folder. All of the functions source code is contained within the src directory. The analysis conducted for Figure 4 was completed on an HPC system. As such, we have included the simulation output files used to construct the final figure within the data folder. Original simulation files for the HPC are included within the "cluster" folder.

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

