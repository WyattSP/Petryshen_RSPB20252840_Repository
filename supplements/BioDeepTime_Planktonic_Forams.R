####################
# NeoTome Analysis #
####################
# Run on Cluster 
# Haar fluctuation analysis of empirical data from BioDeepTime
# Time series of beta diversity from BioDeepTime
library(chronosphere)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(zoo)
library(parallel)
library(stringr)
library(jsonlite)
library(divDyn)

#############################
### Imports and Functions ###
#############################
#Hérbert et al. 2021 Function
# Hérbert et al. 2021 Function. This function can be retrieved from: https://zenodo.org/records/5037581
source("../Haar UNTB Project/Ecological_Modelling/RphlHbrt-rscaling-3c91aa69c10c/R/Haar.R")

# Main loop functions
plantonic_haar_analysis_mc <- function(series_list, scaling_factor = 1000000, ncores = 10){
  series_in <- series_list
  iter_l <- length(series_in$seriesID)
  
  mclapply(1:iter_l, function(i){
    series_name <- series_in$seriesID[i]
    
    # Filter & format
    data_series <- bdt %>%
      filter(seriesID == series_name) %>%
      mutate(tax = word(genus, 1))
    
    occ_bin <- data_series %>%
      group_by(age, tax) %>%
      summarise(
        abundance = sum(abundance, na.rm = TRUE),
        occs = n(),
        .groups = "drop"
      )

    d_in <- data.frame(age = occ_bin$age / scaling_factor, species = occ_bin$tax)
    d_in_haar <- divDyn::divDyn(d_in, tax = "species", bin = "age")
    
    # Haar transforms
    # DivRT
    haarH_RT <- Haar(zoo(d_in_haar$divRT, order.by = as.numeric(d_in_haar$age)),
                     scales=c(1,1/2), abs.scales=FALSE, q=1, overlap=0.25, return.flucs=FALSE)
    
    haarH_Ext <- Haar(zoo(d_in_haar$ext2f3, order.by = as.numeric(d_in_haar$age)),
                      scales=c(1,1/2), abs.scales=FALSE, q=1, overlap=0.25, return.flucs=FALSE)
    
    haarH_Ori <- Haar(zoo(d_in_haar$ori2f3, order.by = as.numeric(d_in_haar$age)),
                      scales=c(1,1/2), abs.scales=FALSE, q=1, overlap=0.25, return.flucs=FALSE)
    
    # DivSIB
    haarH_SIB <- Haar(zoo(d_in_haar$divSIB, order.by = as.numeric(d_in_haar$age)),
                      scales=c(1,1/2), abs.scales=FALSE, q=1, overlap=0.25, return.flucs=FALSE)
    
    # Construct result
    out <- list(
      seriesID = series_name,
      unit = series_in$Unit[i],
      group = series_in$group[i],
      db = series_in$db[i],
      diversity_rt = d_in_haar$divRT,
      diversity_sib = d_in_haar$divSIB,
      age_vector = d_in_haar$age,
      div_rt_fluctuation = haarH_RT$Fluc,
      ext_fluctuation = haarH_Ext$Fluc,
      ori_fluctuation = haarH_Ori$Fluc,
      div_sib_fluctuation = haarH_SIB$Fluc,
      scales = haarH_RT$Scales,
      freq = haarH_RT$Freq,
      Age = haarH_RT$Age,
      lat = data_series$lat[1],
      long = data_series$long[1],
      env = data_series$environment[1]
    )
    
    # Save result to JSON
    save_dir <- "../Plantonic_Results"
    dir.create(save_dir, showWarnings = FALSE)
    save_name <- paste0("neotoma_", series_name, "_results.json")
    save_path_full <- file.path(save_dir, save_name)
    jsonlite::write_json(out, save_path_full, pretty = TRUE)
    
    # Print progress
    message(paste0("Series ", i, " (", series_name, ") saved to JSON"))
    
  }, mc.cores = ncores)
  message("Analysis Complete")
}

#############################
### Load BioDeepTime Data ###
#############################
# download data, verbose=FALSE hides the default chatter 
bdt <- fetch("biodeeptime", verbose=FALSE)

#############################
# Filter and Retrieve Data  #
#############################
# Read RDS of time series diagnostics
series_lengths <- readRDS("../BioDeepTime/series_iter_BDT_FIG.rds")
# Time series with > 100 points
l_series <- series_lengths[series_lengths$nPoints > 10,] # Minimum of 10 time points

#########
#  Run  #
#########
# Filter for only Planktonic Foraminfera
planktonic_forams <- l_series %>% filter(group == "Planktonic foraminifers")

# Run Analysis
plantonic_haar_analysis_mc(planktonic_forams)


