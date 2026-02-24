# Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis

This repository contains the code needed to perform the simulation experiments and analysis described in the following paper: 

Detecting Stabilizing Dynamics in Biased Biodiversity Time Series using Haar Fluctuation Analysis. 

DOI: 10.1098/rspb.2025.2840

Petryshen, W. 1, Pincelli, HM. 1, Vasseur, D. 2

1. Department of Earth and Planetary Sciences, Yale University, New Haven, CT, USA.
2. Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA.

Please make sure all directory locations and save folders are correctly set within each notebook prior to running analysis. The repository also contains a requirements.txt file to ensure all analysis dependencies can be recreated. 

Data needs to be unzip prior to running analysis.

Main correspondence should be directed to Wyatt Petryshen (wyatt.petryshen@yale.edu).

## Abstract 

Characterizing how biodiversity has changed through Earth’s history and uncovering the processes that have driven those changes remain a significant challenge. Haar fluctuation analysis, a recently developed time series method, has been suggested as a powerful new tool to infer macroevolutionary drivers and assess system stability. Yet the ability of this method to identify unique drivers or the timing and dynamics of biodiversity, particularly in biased time series, has not been demonstrated. Here, we assess Haar fluctuation analysis and cross-Haar correlations using process-based ecological simulations that incorporate realistic sampling and depositional biases. We find that simpler (neutral) mechanisms can produce patterns observed in the Phanerozoic record, and that uneven sampling and sedimentary hiatuses can distort scaling relationships, cautioning against mechanistic interpretations. Nonetheless, Haar fluctuation analysis can reliably distinguish stabilizing from non-stabilizing dynamics, even under severe sampling bias, supporting the identification of a long-term equilibrium in Phanerozoic marine biodiversity. Our results suggest that this aspect of Haar fluctuation analysis will be consistent for all cases where the systems return time is longer than the temporal grain of the sample and shorter than the length of the time series, and highlights the value of timescale-based approaches for studying biodiversity dynamics.

## Getting Started

This repository is organized such that each publication figure can be recreated within individual Juypter notebooks. All of the functions source code is additionally contained within the src directory. The analysis conducted for Figure 4 was completed on an HPC system. As such, we have included the simulation output files used to construct the final figure within the data folder. Original simulation files for the HPC are included within the "cluster" folder.

Plots included within the supplementary information are contained at the end of the Fig4_Analysis.ipynb notebook or within the Supplementary Analysis notebooks.

Directory of repository for primary analysis and plotting:
* Figure 2 Analysis and Plot: Fig2_Analysis.ipynb
* Figure 3 Analysis and Plot: Fig3_Analysis.ipynb
* Figure 4 Plot: Fig4_Analysis.ipynb

Additioanl Files:
* Source code: src folder
* Supplementary Analysis:
    * Depositional_Model.ipynb
    * Haar_Algorithm_Comparison.ipynb
    * Turnover_Dynamics.ipynb
    * UNTB_Parameter_Choice.ipynb
    * divDynAnalysis.R
    * Exogenous_Factor_Example.ipynb
    * BioDeepTime_Analysis.ipynb
