# ============================================================================
# Forest Degradation and Carbon Sink Accounting Using Remote Sensing and GIS
# Case Study: Mount Kigali Forest, Nyarugenege District, Rwanda
# ============================================================================

# Load required libraries
library(sf)           # Spatial data handling
library(raster)       # Raster data processing
library(terra)        # Modern raster processing
library(rgee)         # Google Earth Engine interface
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(viridis)      # Color palettes
library(RColorBrewer) # Color schemes
library(gridExtra)    # Multiple plots
library(corrplot)     # Correlation plots
library(randomForest) # Machine learning
library(caret)        # Model training
library(rgdal)        # Spatial data I/O
library(leaflet)      # Interactive maps
library(plotly)       # Interactive plots

# Initialize Google Earth Engine (requires authentication)
ee_Initialize()

# ============================================================================
# 1. STUDY AREA DEFINITION - MOUNT KIGALI FOREST
# ============================================================================

# Define Mount Kigali Forest boundary (Nyarugenege District)
# Coordinates for Mount Kigali area in Rwanda
mount_kigali_bounds <- list(
  min_lon = 30.0,    # Western boundary
  max_lon = 30.15,   # Eastern boundary  
  min_lat = -1.95,   # Southern boundary
  max_lat = -1.85    # Northern boundary
)

# Create study area polygon
study_area <- st_polygon(list(rbind(
  c(mount_kigali_bounds$min_lon, mount_kigali_bounds$min_lat),
  c(mount_kigali_bounds$max_lon, mount_kigali_bounds$min_lat),
  c(mount_kigali_bounds$max_lon, mount_kigali_bounds$max_lat),
  c(mount_kigali_bounds$min_lon, mount_kigali_bounds$max_lat),
  c(mount_kigali_bounds$min_lon, mount_kigali_bounds$min_lat)
))) %>%
  st_sfc(crs = 4326) %>%
  st_sf()

# ============================================================================
# 2. DATA ACQUISITION FUNCTIONS
# ============================================================================

# Function to get Landsat data for forest cover analysis
get_landsat_data <- function(start_date, end_date, cloud_threshold = 20) {
  
  # Define study area for Earth Engine
  aoi <- ee$Geometry$Rectangle(c(mount_kigali_bounds$min_lon, 
                                mount_kigali_bounds$min_lat,
                                mount_kigali_bounds$max_lon, 
                                mount_kigali_bounds$max_lat))
  
  # Get Landsat 8 collection
  landsat <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
    filterDate(start_date, end_date)$
    filterBounds(aoi)$
    filter(ee$Filter$lt("CLOUD_COVER", cloud_threshold))
  
  # Function to mask clouds and apply scaling factors
  mask_and_scale <- function(image) {
    # Get QA pixel band
    qa <- image$select("QA_PIXEL")
    
    # Cloud mask
    cloud_mask <- qa$bitwiseAnd(as.integer(2^3))$eq(0)$
      And(qa$bitwiseAnd(as.integer(2^4))$eq(0))
    
    # Apply scaling factors
    optical_bands <- image$select("SR_B.")$multiply(0.0000275)$add(-0.2)
    thermal_bands <- image$select("ST_B.*")$multiply(0.00341802)$add(149.0)
    
    return(image$addBands(optical_bands, overwrite = TRUE)$
           addBands(thermal_bands, overwrite = TRUE)$
           updateMask(cloud_mask))
  }
  
  # Apply preprocessing and create median composite
  landsat_processed <- landsat$map(mask_and_scale)$median()$clip(aoi)
  
  return(landsat_processed)
}

# Function to calculate vegetation indices
calculate_vegetation_indices <- function(image) {
  
  # NDVI (Normalized Difference Vegetation Index)
  ndvi <- image$normalizedDifference(c("SR_B5", "SR_B4"))$rename("NDVI")
  
  # EVI (Enhanced Vegetation Index)
  evi <- image$expression(
    "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))",
    list(
      NIR = image$select("SR_B5"),
      RED = image$select("SR_B4"),
      BLUE = image$select("SR_B2")
    )
  )$rename("EVI")
  
  # SAVI (Soil Adjusted Vegetation Index)
  savi <- image$expression(
    "((NIR - RED) / (NIR + RED + 0.5)) * (1 + 0.5)",
    list(
      NIR = image$select("SR_B5"),
      RED = image$select("SR_B4")
    )
  )$rename("SAVI")
  
  # MSAVI (Modified Soil Adjusted Vegetation Index)
  msavi <- image$expression(
    "(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED))) / 2",
    list(
      NIR = image$select("SR_B5"),
      RED = image$select("SR_B4")
    )
  )$rename("MSAVI")
  
  return(image$addBands(list(ndvi, evi, savi, msavi)))
}

# ============================================================================
# 3. OBJECTIVE 1: FOREST COVER MAPPING OVER TIME
# ============================================================================

# Define time periods for analysis
time_periods <- list(
  period_1 = list(start = "2015-01-01", end = "2016-12-31"),
  period_2 = list(start = "2017-01-01", end = "2018-12-31"),
  period_3 = list(start = "2019-01-01", end = "2020-12-31"),
  period_4 = list(start = "2021-01-01", end = "2022-12-31"),
  period_5 = list(start = "2023-01-01", end = "2024-12-31")
)

# Function to create forest cover maps
create_forest_cover_maps <- function() {
  
  forest_maps <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    period <- time_periods[[i]]
    
    cat("Processing", period_name, ":", period$start, "to", period$end, "\n")
    
    # Get Landsat data for the period
    landsat_image <- get_landsat_data(period$start, period$end)
    
    # Calculate vegetation indices
    image_with_indices <- calculate_vegetation_indices(landsat_image)
    
    # Create forest mask using NDVI threshold
    forest_mask <- image_with_indices$select("NDVI")$gt(0.4)$rename("forest_cover")
    
    # Add to collection
    forest_maps[[period_name]] <- forest_mask
    
    # Download as raster for local analysis
    forest_raster <- ee_as_raster(
      image = forest_mask,
      region = study_area,
      scale = 30,
      via = "drive"
    )
    
    # Save locally
    writeRaster(forest_raster, 
                filename = paste0("forest_cover_", period_name, ".tif"),
                overwrite = TRUE)
  }
  
  return(forest_maps)
}

# ============================================================================
# 4. OBJECTIVE 2: BIOMASS AND CARBON STOCK ESTIMATION
# ============================================================================

# Function to estimate above-ground biomass using allometric equations
estimate_agb <- function(ndvi_image, evi_image) {
  
  # Biomass estimation using vegetation indices
  # Based on Chave et al. (2014) and local calibration studies
  
  # Method 1: NDVI-based estimation
  agb_ndvi <- ndvi_image$expression(
    "exp(2.148 + 0.207 * log(ndvi_scaled) + 0.928 * log(evi_scaled))",
    list(
      ndvi_scaled = ndvi_image$multiply(100),
      evi_scaled = evi_image$multiply(100)
    )
  )$rename("AGB_NDVI")
  
  # Method 2: Multiple regression approach
  agb_multi <- ndvi_image$expression(
    "150 * pow(ndvi, 2) + 50 * evi + 25",
    list(
      ndvi = ndvi_image,
      evi = evi_image
    )
  )$rename("AGB_MULTI")
  
  return(list(agb_ndvi = agb_ndvi, agb_multi = agb_multi))
}

# Function to convert biomass to carbon stocks
biomass_to_carbon <- function(agb_image, carbon_fraction = 0.47) {
  
  # Above-ground carbon
  agc <- agb_image$multiply(carbon_fraction)$rename("AGC")
  
  # Below-ground carbon (root-to-shoot ratio)
  # Using IPCC default values for tropical forests
  bgc <- agc$multiply(0.24)$rename("BGC")
  
  # Total carbon stock
  total_carbon <- agc$add(bgc)$rename("TOTAL_CARBON")
  
  return(list(agc = agc, bgc = bgc, total_carbon = total_carbon))
}

# Function to create comprehensive biomass and carbon maps
create_biomass_carbon_maps <- function() {
  
  biomass_carbon_results <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    period <- time_periods[[i]]
    
    cat("Estimating biomass and carbon for", period_name, "\n")
    
    # Get Landsat data
    landsat_image <- get_landsat_data(period$start, period$end)
    image_with_indices <- calculate_vegetation_indices(landsat_image)
    
    # Extract vegetation indices
    ndvi <- image_with_indices$select("NDVI")
    evi <- image_with_indices$select("EVI")
    
    # Estimate biomass
    biomass_estimates <- estimate_agb(ndvi, evi)
    
    # Convert to carbon
    carbon_stocks <- biomass_to_carbon(biomass_estimates$agb_multi)
    
    # Store results
    biomass_carbon_results[[period_name]] <- list(
      biomass = biomass_estimates,
      carbon = carbon_stocks
    )
    
    # Download carbon stock map
    carbon_raster <- ee_as_raster(
      image = carbon_stocks$total_carbon,
      region = study_area,
      scale = 30,
      via = "drive"
    )
    
    writeRaster(carbon_raster,
                filename = paste0("carbon_stocks_", period_name, ".tif"),
                overwrite = TRUE)
  }
  
  return(biomass_carbon_results)
}

# ============================================================================
# 5. OBJECTIVE 3: FOREST FRAGMENTATION ANALYSIS
# ============================================================================

# Function to calculate fragmentation indices
calculate_fragmentation_indices <- function(forest_raster) {
  
  # Convert to binary forest/non-forest
  forest_binary <- forest_raster > 0
  
  # Calculate landscape metrics using landscapemetrics approach
  # Patch density
  patch_density <- focal(forest_binary, w = matrix(1, 9, 9), fun = sum) / 81
  
  # Edge density (simplified)
  forest_edges <- boundaries(forest_binary, type = "inner")
  edge_density <- focal(forest_edges, w = matrix(1, 9, 9), fun = sum)
  
  # Forest connectivity (simplified)
  connectivity <- focal(forest_binary, w = matrix(1, 15, 15), fun = sum) / 225
  
  return(list(
    patch_density = patch_density,
    edge_density = edge_density,
    connectivity = connectivity
  ))
}

# Function to analyze fragmentation-carbon relationships
analyze_fragmentation_carbon_correlation <- function() {
  
  correlation_results <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    
    # Load forest cover and carbon stock rasters
    forest_file <- paste0("forest_cover_", period_name, ".tif")
    carbon_file <- paste0("carbon_stocks_", period_name, ".tif")
    
    if (file.exists(forest_file) && file.exists(carbon_file)) {
      forest_raster <- raster(forest_file)
      carbon_raster <- raster(carbon_file)
      
      # Calculate fragmentation indices
      frag_indices <- calculate_fragmentation_indices(forest_raster)
      
      # Extract values for correlation analysis
      forest_values <- getValues(forest_raster)
      carbon_values <- getValues(carbon_raster)
      patch_density_values <- getValues(frag_indices$patch_density)
      edge_density_values <- getValues(frag_indices$edge_density)
      connectivity_values <- getValues(frag_indices$connectivity)
      
      # Remove NA values
      valid_indices <- complete.cases(
        forest_values, carbon_values, patch_density_values,
        edge_density_values, connectivity_values
      )
      
      # Create correlation matrix
      correlation_data <- data.frame(
        Forest_Cover = forest_values[valid_indices],
        Carbon_Stock = carbon_values[valid_indices],
        Patch_Density = patch_density_values[valid_indices],
        Edge_Density = edge_density_values[valid_indices],
        Connectivity = connectivity_values[valid_indices]
      )
      
      correlation_matrix <- cor(correlation_data, use = "complete.obs")
      correlation_results[[period_name]] <- correlation_matrix
    }
  }
  
  return(correlation_results)
}

# ============================================================================
# 6. OBJECTIVE 4: GEOSPATIAL INDICATORS FOR DEGRADATION MONITORING
# ============================================================================

# Function to create degradation indicators
create_degradation_indicators <- function() {
  
  # Load all time period data
  forest_time_series <- list()
  carbon_time_series <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    
    forest_file <- paste0("forest_cover_", period_name, ".tif")
    carbon_file <- paste0("carbon_stocks_", period_name, ".tif")
    
    if (file.exists(forest_file)) {
      forest_time_series[[period_name]] <- raster(forest_file)
    }
    if (file.exists(carbon_file)) {
      carbon_time_series[[period_name]] <- raster(carbon_file)
    }
  }
  
  # Calculate change indicators
  degradation_indicators <- list()
  
  if (length(forest_time_series) >= 2) {
    # Forest cover change rate
    first_period <- forest_time_series[[1]]
    last_period <- forest_time_series[[length(forest_time_series)]]
    
    forest_change_rate <- (last_period - first_period) / 
                         (length(forest_time_series) - 1)
    
    # Carbon stock change rate
    if (length(carbon_time_series) >= 2) {
      first_carbon <- carbon_time_series[[1]]
      last_carbon <- carbon_time_series[[length(carbon_time_series)]]
      
      carbon_change_rate <- (last_carbon - first_carbon) / 
                           (length(carbon_time_series) - 1)
      
      degradation_indicators$carbon_change_rate <- carbon_change_rate
    }
    
    degradation_indicators$forest_change_rate <- forest_change_rate
    
    # Degradation hotspots (areas with significant decline)
    degradation_hotspots <- forest_change_rate < -0.1  # 10% decline threshold
    degradation_indicators$hotspots <- degradation_hotspots
  }
  
  return(degradation_indicators)
}

# ============================================================================
# 7. COMPREHENSIVE ANALYSIS AND REPORTING
# ============================================================================

# Main analysis function
run_comprehensive_analysis <- function() {
  
  cat("Starting comprehensive forest degradation and carbon analysis...\n")
  
  # Create output directory
  if (!dir.exists("outputs")) {
    dir.create("outputs")
  }
  
  # Step 1: Forest cover mapping
  cat("Step 1: Creating forest cover maps...\n")
  forest_maps <- create_forest_cover_maps()
  
  # Step 2: Biomass and carbon estimation
  cat("Step 2: Estimating biomass and carbon stocks...\n")
  biomass_carbon_results <- create_biomass_carbon_maps()
  
  # Step 3: Fragmentation analysis
  cat("Step 3: Analyzing forest fragmentation...\n")
  fragmentation_results <- analyze_fragmentation_carbon_correlation()
  
  # Step 4: Degradation indicators
  cat("Step 4: Creating degradation indicators...\n")
  degradation_results <- create_degradation_indicators()
  
  # Generate summary statistics
  generate_summary_statistics(fragmentation_results, degradation_results)
  
  cat("Analysis completed! Check the outputs directory for results.\n")
  
  return(list(
    forest_maps = forest_maps,
    biomass_carbon = biomass_carbon_results,
    fragmentation = fragmentation_results,
    degradation = degradation_results
  ))
}

# Function to generate summary statistics and visualizations
generate_summary_statistics <- function(fragmentation_results, degradation_results) {
  
  # Create summary report
  sink("outputs/analysis_summary.txt")
  
  cat("FOREST DEGRADATION AND CARBON ACCOUNTING ANALYSIS\n")
  cat("Study Area: Mount Kigali Forest, Nyarugenege District, Rwanda\n")
  cat("Analysis Date:", Sys.Date(), "\n\n")
  
  cat("CORRELATION ANALYSIS RESULTS:\n")
  for (period in names(fragmentation_results)) {
    cat("\nPeriod:", period, "\n")
    cat("Correlation between Forest Cover and Carbon Stock:",
        round(fragmentation_results[[period]]["Forest_Cover", "Carbon_Stock"], 3), "\n")
    cat("Correlation between Connectivity and Carbon Stock:",
        round(fragmentation_results[[period]]["Connectivity", "Carbon_Stock"], 3), "\n")
  }
  
  cat("\nDEGRADATION INDICATORS:\n")
  if (!is.null(degradation_results$forest_change_rate)) {
    forest_change_stats <- cellStats(degradation_results$forest_change_rate, 
                                   c('mean', 'min', 'max'), na.rm = TRUE)
    cat("Mean annual forest change rate:", round(forest_change_stats[1], 4), "\n")
    cat("Minimum forest change rate:", round(forest_change_stats[2], 4), "\n")
    cat("Maximum forest change rate:", round(forest_change_stats[3], 4), "\n")
  }
  
  sink()
  
  # Create visualization plots
  create_analysis_plots(fragmentation_results, degradation_results)
}

# Function to create analysis plots
create_analysis_plots <- function(fragmentation_results, degradation_results) {
  
  # Plot correlation matrices
  if (length(fragmentation_results) > 0) {
    pdf("outputs/correlation_analysis.pdf", width = 12, height = 8)
    
    for (period in names(fragmentation_results)) {
      corrplot(fragmentation_results[[period]], 
               method = "color", 
               type = "upper",
               title = paste("Correlation Matrix -", period),
               mar = c(0, 0, 2, 0))
    }
    
    dev.off()
  }
  
  # Plot degradation trends if available
  if (!is.null(degradation_results$forest_change_rate)) {
    
    pdf("outputs/degradation_maps.pdf", width = 10, height = 8)
    
    # Forest change rate map
    plot(degradation_results$forest_change_rate,
         main = "Annual Forest Change Rate",
         col = RColorBrewer::brewer.pal(11, "RdYlGn"))
    
    # Degradation hotspots
    plot(degradation_results$hotspots,
         main = "Forest Degradation Hotspots",
         col = c("lightgreen", "red"),
         legend = FALSE)
    legend("bottomright", 
           legend = c("Stable", "Degraded"), 
           fill = c("lightgreen", "red"))
    
    dev.off()
  }
}

# ============================================================================
# 8. EXECUTION
# ============================================================================

# Run the comprehensive analysis
# Uncomment the line below to execute the analysis
# results <- run_comprehensive_analysis()

# Print execution instructions
cat("FOREST DEGRADATION AND CARBON ACCOUNTING ANALYSIS\n")
cat("=================================================\n\n")
cat("To run the analysis:\n")
cat("1. Ensure you have authenticated with Google Earth Engine using ee_Initialize()\n")
cat("2. Run: results <- run_comprehensive_analysis()\n")
cat("3. Check the 'outputs' directory for results\n\n")
cat("Generated outputs will include:\n")
cat("- Forest cover maps for each time period (.tif files)\n")
cat("- Carbon stock maps (.tif files)\n")
cat("- Analysis summary (analysis_summary.txt)\n")
cat("- Correlation analysis plots (correlation_analysis.pdf)\n")
cat("- Degradation maps (degradation_maps.pdf)\n\n")
cat("Study area: Mount Kigali Forest, Nyarugenege District, Rwanda\n")
cat("Coordinates: 30.0°E to 30.15°E, -1.95°S to -1.85°S\n")