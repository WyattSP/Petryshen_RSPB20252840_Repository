# Added as a revision
# R Script to extract BioDeepTime data and conduct Haar Fluctuation Analysis

# Haar fluctuation analysis of empirical data from BioDeepTime
# Time series of beta diversity from BioDeepTime
library(chronosphere)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(zoo)
library(viridis)
library(maps)
library(viridis) 
library(ncdf4)
library(jsonlite)

#############################
### Import Haar Functions ###
#############################
# Petryshen Haar Function
#source("/Users/wyattpetryshen/Documents/UNTB_HAAR/code/haar_wp.R")
# Haar Fluc from Arnischidt
source("../Haar_fluc.R")
# Hérbert et al. 2021 Function. This function can be retrieved from: https://zenodo.org/records/5037581
source("../RphlHbrt-rscaling-3c91aa69c10c/R/Haar.R")

#############################
### Functions Definitions ###
#############################
# Plot Function
plot_haar_all <- function(haar_output_list){
  # Mins
  min_ages <- sapply(haar_output_list, function(x) min(x$Age, na.rm = TRUE))
  age_colors <- viridis(length(min_ages))[rank(min_ages)]  
  
  # Initialize Plot
  plot(NA, NA, xlim = c(-2, 8), ylim = c(-1.5, 1.5), 
       xlab = "log10(Scale)", ylab = "log10(Fluctuation)", 
       main = "", type = "n")
  
  # Plot all
  for (i in seq_along(haar_output_list)) {
    # Extract series metadata and haar output
    item <- haar_output_list[[i]]
    series_id <- item$seriesID
    
    # Compute min age
    min_age <- min(item$Age, na.rm = TRUE)
    
    # Create plot data frame
    df <- data.frame(
      fluctuation = item$fluctuation,
      scales = item$scales
    )
    lines(log10(df$scales), log10(df$fluctuation), col = age_colors[i], lwd = 2, main = item$seriesID)
  }
}

# Haar analysis
haar_analysis <- function(series_list){
  series_in <- series_list
  haar_results_list <- list()
  for(i in 1:length(series_in$seriesID)){
    # Series ID
    series_name = series_in$seriesID[i]
    # Filter Data
    data_in <- bdt %>% filter(seriesID == series_name)
    # Order data
    result <- data_in %>%
      group_by(age) %>%
      summarise(unique_species_count = n_distinct(analyzedTaxon))
    # Convert to Zoo object
    zooSeries <- zoo(result$unique_species_count, order.by = as.numeric(result$age))
    # Haar Analysis
    haarH <- Haar(zooSeries, scales=NULL, abs.scales=FALSE, q=1, overlap=0.5, return.flucs=FALSE)
    # Extract Haar outputs and convert to data frame
    haar_df <- list(
      seriesID = series_name,
      unit = series_in$Unit[i],
      group = series_in$group[i],
      db = series_in$db[i],
      fluctuation = haarH$Fluc,
      scales = haarH$Scales,
      freq = haarH$Freq,
      SDHaar = haarH$SDHaarFluc,
      Age = haarH$Age,
      lat = data_in$lat[1],
      long = data_in$long[1],
      env = data_in$environment[1]
    )
    # Append to results list
    haar_results_list[[i]] <- haar_df
  }
  return(haar_results_list)
}

# Slope Estimator
slope_estimator <- function(haar_results_list, range_min, range_max){
  # We are assuming at least 2 break points within each Haar fluctuation
  breakpoints_haar = list()
  # Range to estimate slope
  r_min <- range_min
  r_max <- range_max
  
  # Loop to estimate
  for(i in 1:length(haar_results_list)){
    # Get data
    item <- haar_results_list[[i]]
    x = log10(item$scales)
    y = log10(item$fluc) 
    
    # Filter based on range
    mask <- x >= r_min & x <= r_max
    x <- x[mask]
    y <- y[mask]
    
    # Check for NaN or NA in y (fluctuations)
    if (any(is.nan(y)) || any(is.na(y))) {
      message(sprintf("Skipping iteration %d due to NaN/NA in fluctuations", i))
      breakpoints_haar[[i]] <- NA
      next  # Skip to next iteration
    }
    
    # Fit naive lm
    lm.out <- lm(y ~ x)
    #summary(lm.out)
    my.coef <- coef(lm.out)
    
    # Using segmented and breakpoint analysis
    my.seg <- segmented(lm.out,
                        seg.Z = ~ x,
                        psi = list(x = median(x, na.rm = TRUE))) # Weak point is estimate for psi
    
    # Extract relevant output
    breakpoint_out <- list(
      seriesID = item$seriesID,
      unit = item$unit,
      group = item$group,
      lat = item$lat[1],
      long = item$long[1],
      long = item$env[1],
      breaks_out = my.seg$psi,
      slopes_out = slope(my.seg),
      segments = my.seg,
      max_scale = max(x),
      x = x,
      y = y
    )
    # Append to results list
    breakpoints_haar[[i]] <- breakpoint_out
  }
  return(breakpoints_haar)
}

# Extract Values #
extract_slope_info <- function(slopes_list) {
  list(
    s1 = unlist(lapply(slopes_list, function(x) x$slopes_out$x[1])),
    s2 = unlist(lapply(slopes_list, function(x) x$slopes_out$x[2])),
    bp1 = unlist(lapply(slopes_list, function(x) x$breaks_out[2])),
    max_x = unlist(lapply(slopes_list, function(x) max(x$x))),
    median_x = unlist(lapply(slopes_list, function(x) median(x$x))),
    lat = unlist(lapply(slopes_list, function(x) x$lat)),
    long = unlist(lapply(slopes_list, function(x) x$long)),
    env = unlist(lapply(slopes_list, function(x) x$env)),
    min_age = unlist(lapply(slopes_list, function(x) 10^x$x[1]))
  )
}

# Plot Point Map #
plot_point_map <- function(long, lat, value){
  world_map <- map_data("world")
  
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "gray60") +
    geom_point(aes(x = (long), y = (lat), color = (value)), size = 5) +
    scale_color_viridis_c(option = "plasma", 
                          name = "Signal Intensity") +
    theme_minimal() +
    labs(title = "", x = "Longitude", y = "Latitude") +
    theme(legend.position = "right")
}

# Grid Function
calculate_grid_values <- function(lon, lat, var, grid_long, grid_lat) {
  mesh_grid <- expand.grid(lon = grid_long, lat = grid_lat)
  mesh_grid$mean_val <- NA
  mesh_grid$sd_val <- NA
  
  lon_res <- diff(grid_long)[1]
  lat_res <- diff(grid_lat)[1]
  
  for (i in 1:nrow(mesh_grid)) {
    cell_long <- mesh_grid$lon[i]
    cell_lat <- mesh_grid$lat[i]
    
    mask <- lon >= (cell_long - lon_res/2) & lon < (cell_long + lon_res/2) &
      lat >= (cell_lat - lat_res/2) & lat < (cell_lat + lat_res/2)
    
    if (any(mask)) {
      mesh_grid$mean_val[i] <- mean(var[mask], na.rm = TRUE)
      mesh_grid$sd_val[i] <- sd(var[mask], na.rm = TRUE)
    }
  }
  
  return(mesh_grid)
}

# Export Function
export_as_netcdf <- function(grid_data_values, grid_long, grid_lat, var_name, save_path){
  # Reshape the result into a matrix (lon x lat), respecting dimension order in NetCDF
  mean_matrix <- matrix(grid_data_values, 
                        nrow = length(grid_long), 
                        ncol = length(grid_lat), 
                        byrow = FALSE)
  # Define NetCDF dimensions (make sure names match!)
  dim_lon <- ncdim_def("longitude", units = "degrees_east", vals = grid_long)
  dim_lat <- ncdim_def("latitude", units = "degrees_north", vals = grid_lat)
  # Define variable
  var_def <- ncvar_def(var_name, units = "unitless", dim = list(dim_lon, dim_lat), 
                       missval = NA, prec = "double")
  # Save as NC
  nc_out <- nc_create(save_path, vars = var_def)
  ncvar_put(nc_out, var_def, mean_matrix)
  nc_close(nc_out)
  print("Save Complete")
}

#############################
### Load BioDeepTime Data ###
#############################
# download data, verbose=FALSE hides the default chatter 
bdt <- fetch("biodeeptime", verbose=FALSE)

#############################
### Load BioDeepTime Data ###
#############################
# download data, verbose=FALSE hides the default chatter 
bdt <- fetch("biodeeptime", verbose=FALSE)

# Modern
modern <- bdt %>%
  filter(db == "BioTIME")

#############################
# Filter and Retrieve Data  #
#############################
# Read RDS of time series diagnostics
series_lengths <- readRDS("../BioDeepTime/series_iter_BDT_FIG.rds")

# Time series with > 100 points
l_series <- series_lengths[series_lengths$nPoints > 50,] # Minimum of 50 time points

# Filter based on time units
# Units ma
ma_tseries <- l_series[grepl("ma", l_series$Unit), ]
# Units ka
ka_tseries <- l_series[grepl("ka", l_series$Unit), ]
# Units Years 
yr_tseries <- l_series[grepl("years", l_series$Unit), ]
# Units days
day_tseries <- l_series[grepl("rday", l_series$Unit), ]

####################################
### Analysis of MA Filtered Data ###
####################################
# Haar Analysis #
ma_haar_results_list <- haar_analysis(ma_tseries)
# Diagnostic Plot #
plot_haar_all(ma_haar_results_list)
# Slope Estimation #
ma_slopes <- slope_estimator(ma_haar_results_list, log10(2), log10(100e6))
# Remove NAN #
ma_slopes_clean <- Filter(function(x) {!is.na(x[[1]])}, ma_slopes)
# Extract Values #
ma_slope_info <- extract_slope_info(ma_slopes_clean)
# Plot #
plot_point_map(ma_slope_info$long, ma_slope_info$lat, ma_slope_info$s1)
# Grid to HadCM3 Resolution #
grid_long <- seq(from = 0, to = 356.25, by = 3.75)
grid_lat <- seq(from = -90, to = 90, by = 2.5)
long <- ma_slope_info$long
long[long < 0] <- long[long < 0] + 360
ma_grid_data <- calculate_grid_values(long, ma_slope_info$lat, ma_slope_info$s1, grid_long, grid_lat)
# Export as NetCDF4 #
save_dir <- "../NC_Haar_Beta_Files"
save_path_full <- file.path(save_dir, "ma_beta.nc")
export_as_netcdf(ma_grid_data$mean_val, grid_long, grid_lat, "beta", save_path_full)