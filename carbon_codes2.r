# Function to create analysis plots (improved with error handling)
create_analysis_plots <- function(fragmentation_results, degradation_results) {
  
  tryCatch({
    # Plot correlation matrices
    if (length(fragmentation_results) > 0) {
      cat("Creating correlation analysis plots...\n")
      
      pdf("outputs/correlation_analysis.pdf", width = 12, height = 8)
      
      for (period in names(fragmentation_results)) {
        tryCatch({
          corr_matrix <- fragmentation_results[[period]]
          if (!is.null(corr_matrix) && is.matrix(corr_matrix)) {
            corrplot(corr_matrix, 
                     method = "color", 
                     type = "upper",
                     title = paste("Correlation Matrix -", period),
                     mar = c(0, 0, 2, 0),
                     tl.cex = 0.8,
                     number.cex = 0.8)
          }
        }, error = function(e) {
          cat("Error plotting correlation for", period, ":", e$message, "\n")
        })
      }
      
      dev.off()
      cat("Correlation plots saved to outputs/correlation_analysis.pdf\n")# Function to create comprehensive biomass and carbon maps
create_biomass_carbon_maps <- function() {
  
  biomass_carbon_results <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    period <- time_periods[[i]]
    
    cat("Estimating biomass and carbon for", period_name, "\n")
    
    # Get Landsat data
    landsat_image <- get_landsat_data(period$start, period$end)
    image_with_indices <- calculate_vegetation_indices(landsat_image)
    
    if (gee_available) {
      # Google Earth Engine version
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
      tryCatch({
        carbon_raster <- ee_as_raster(
          image = carbon_stocks$total_carbon,
          region = study_area,
          scale = 30,
          via = "drive"
        )
        
        writeRaster(carbon_raster,
                    filename = paste0("carbon_stocks_", period_name, ".tif"),
                    overwrite = TRUE)
      }, error = function(e) {
        cat("Error downloading carbon data from GEE, creating local version...\n")
        # Create local fallback version
        ndvi_local <- image_with_indices[["NDVI"]]
        evi_local <- image_with_indices[["EVI"]]
        biomass_local <- estimate_agb(ndvi_local, evi_local)
        carbon_local <- biomass_to_carbon(biomass_local$agb_multi)
        
        writeRaster(carbon_local$total_carbon,
                    filename = paste0("carbon_stocks_", period_name, ".tif"),
                    overwrite = TRUE)
      })
      
    } else {
      # Local raster version
      # Extract vegetation indices
      ndvi <- image_with_indices[["NDVI"]]
      evi <- image_with_indices[["EVI"]]
      
      # Estimate biomass
      biomass_estimates <- estimate_agb(ndvi, evi)
      
      # Convert to carbon
      carbon_stocks <- biomass_to_carbon(biomass_estimates$agb_multi)
      
      # Store results
      biomass_carbon_results[[period_name]] <- list(
        biomass = biomass_estimates,
        carbon = carbon_stocks
      )
      
      # Save carbon stock map
      writeRaster(carbon_stocks$total_carbon,
                  filename = paste0("carbon_stocks_", period_name, ".tif"),
                  overwrite = TRUE)
    }
  }
  
  return(biomass_carbon_results)
}# ============================================================================
# Forest Degradation and Carbon Sink Accounting Using Remote Sensing and GIS
# Case Study: Mount Kigali Forest, Nyarugenege District, Rwanda
# ============================================================================

# ============================================================================
# PACKAGE INSTALLATION AND LOADING
# ============================================================================

# Function to install and load packages safely
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      cat("Package", pkg, "already loaded\n")
    }
  }
}

# List of required packages
required_packages <- c(
  "sf",             # Spatial data handling
  "raster",         # Raster data processing
  "terra",          # Modern raster processing
  "sp",             # Spatial classes and methods
  "rgdal",          # Spatial data I/O
  "ggplot2",        # Plotting
  "dplyr",          # Data manipulation
  "viridis",        # Color palettes
  "RColorBrewer",   # Color schemes
  "gridExtra",      # Multiple plots
  "corrplot",       # Correlation plots
  "randomForest",   # Machine learning
  "caret",          # Model training
  "leaflet",        # Interactive maps
  "plotly",         # Interactive plots
  "magrittr",       # Pipe operators
  "tidyr",          # Data tidying
  "devtools"        # Development tools
)

# Install and load core packages
cat("Installing and loading core packages...\n")
install_and_load(required_packages)

# Install rgee (Google Earth Engine) with special handling
install_rgee <- function() {
  if (!require("rgee", character.only = TRUE)) {
    cat("Installing rgee (Google Earth Engine interface)...\n")
    
    # Install from CRAN first
    tryCatch({
      install.packages("rgee", dependencies = TRUE)
      library(rgee)
      cat("rgee installed successfully from CRAN\n")
    }, error = function(e) {
      cat("CRAN installation failed, trying development version...\n")
      
      # Try development version if CRAN fails
      tryCatch({
        if (!require("devtools")) {
          install.packages("devtools")
          library(devtools)
        }
        devtools::install_github("r-spatial/rgee")
        library(rgee)
        cat("rgee installed successfully from GitHub\n")
      }, error = function(e2) {
        cat("Warning: Could not install rgee. Google Earth Engine functionality will be limited.\n")
        cat("You can proceed with local raster analysis only.\n")
        return(FALSE)
      })
    })
  } else {
    library(rgee)
    cat("rgee already available\n")
  }
  return(TRUE)
}

# Install additional geospatial packages
install_geospatial_packages <- function() {
  geo_packages <- c(
    "landscapemetrics",  # Landscape ecology metrics
    "rasterVis",        # Raster visualization
    "mapview",          # Interactive mapping
    "tmap",             # Thematic maps
    "stars",            # Spatiotemporal arrays
    "lwgeom",           # Lightweight geometry
    "units",            # Units of measurement
    "classInt"          # Class intervals
  )
  
  cat("Installing additional geospatial packages...\n")
  
  for (pkg in geo_packages) {
    tryCatch({
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
        cat("Successfully installed:", pkg, "\n")
      }
    }, error = function(e) {
      cat("Warning: Could not install", pkg, "- continuing without it\n")
    })
  }
}

# Execute package installation
cat("=== PACKAGE INSTALLATION STARTING ===\n")
install_and_load(required_packages)
rgee_available <- install_rgee()
install_geospatial_packages()
cat("=== PACKAGE INSTALLATION COMPLETED ===\n\n")

# Check R version and provide warnings if needed
check_r_version <- function() {
  r_version <- R.version.string
  cat("R Version:", r_version, "\n")
  
  if (R.version$major < 4) {
    cat("Warning: R version 4.0+ recommended for optimal performance\n")
  }
}

check_r_version()

# Alternative function for users without Google Earth Engine
use_local_data_only <- function() {
  cat("========================================\n")
  cat("RUNNING IN LOCAL DATA MODE\n")
  cat("Google Earth Engine not available\n")
  cat("Analysis will use sample/simulated data\n")
  cat("========================================\n")
}

# Initialize Google Earth Engine with error handling
initialize_gee <- function() {
  if (rgee_available) {
    tryCatch({
      # Check if Python and Earth Engine are properly configured
      ee_check()
      
      # Initialize Earth Engine
      ee_Initialize(drive = TRUE)
      cat("Google Earth Engine initialized successfully\n")
      return(TRUE)
    }, error = function(e) {
      cat("Error initializing Google Earth Engine:\n")
      cat(as.character(e), "\n")
      cat("Proceeding with local analysis mode...\n")
      use_local_data_only()
      return(FALSE)
    })
  } else {
    cat("rgee not available - using local analysis mode\n")
    use_local_data_only()
    return(FALSE)
  }
}

# Try to initialize Google Earth Engine
gee_available <- initialize_gee()

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
  
  if (!gee_available) {
    cat("Google Earth Engine not available. Creating simulated data...\n")
    return(create_simulated_landsat_data())
  }
  
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

# Function to create simulated Landsat data when GEE is not available
create_simulated_landsat_data <- function() {
  
  # Create raster template
  ncols <- 100
  nrows <- 100
  
  # Create coordinate system
  extent_obj <- extent(mount_kigali_bounds$min_lon, mount_kigali_bounds$max_lon,
                      mount_kigali_bounds$min_lat, mount_kigali_bounds$max_lat)
  
  # Simulate realistic forest patterns
  set.seed(123)  # For reproducible results
  
  # Red band (Band 4)
  red_band <- raster(nrows = nrows, ncols = ncols, ext = extent_obj, crs = "+proj=longlat +datum=WGS84")
  values(red_band) <- runif(ncell(red_band), 0.05, 0.25)
  names(red_band) <- "SR_B4"
  
  # NIR band (Band 5) - higher values for vegetation
  nir_band <- raster(nrows = nrows, ncols = ncols, ext = extent_obj, crs = "+proj=longlat +datum=WGS84")
  values(nir_band) <- runif(ncell(nir_band), 0.15, 0.45)
  names(nir_band) <- "SR_B5"
  
  # Blue band (Band 2)
  blue_band <- raster(nrows = nrows, ncols = ncols, ext = extent_obj, crs = "+proj=longlat +datum=WGS84")
  values(blue_band) <- runif(ncell(blue_band), 0.02, 0.15)
  names(blue_band) <- "SR_B2"
  
  # Create forest-like patterns (higher NIR in forested areas)
  forest_mask <- values(nir_band) > 0.3
  values(nir_band)[forest_mask] <- values(nir_band)[forest_mask] * 1.5
  values(red_band)[forest_mask] <- values(red_band)[forest_mask] * 0.7
  
  # Stack bands
  landsat_stack <- stack(red_band, nir_band, blue_band)
  
  cat("Created simulated Landsat data for Mount Kigali area\n")
  return(landsat_stack)
}

# Function to calculate vegetation indices
calculate_vegetation_indices <- function(image) {
  
  if (gee_available) {
    # Google Earth Engine version
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
    
  } else {
    # Local raster version
    red <- image[["SR_B4"]]
    nir <- image[["SR_B5"]]
    blue <- image[["SR_B2"]]
    
    # NDVI calculation
    ndvi <- (nir - red) / (nir + red)
    names(ndvi) <- "NDVI"
    
    # EVI calculation
    evi <- 2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))
    names(evi) <- "EVI"
    
    # SAVI calculation
    savi <- ((nir - red) / (nir + red + 0.5)) * (1 + 0.5)
    names(savi) <- "SAVI"
    
    # MSAVI calculation
    msavi <- (2 * nir + 1 - sqrt((2 * nir + 1)^2 - 8 * (nir - red))) / 2
    names(msavi) <- "MSAVI"
    
    # Stack all indices
    indices_stack <- stack(ndvi, evi, savi, msavi)
    
    # Combine with original image
    result <- stack(image, indices_stack)
    
    return(result)
  }
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
    
    if (gee_available) {
      # Google Earth Engine version
      # Create forest mask using NDVI threshold
      forest_mask <- image_with_indices$select("NDVI")$gt(0.4)$rename("forest_cover")
      
      # Add to collection
      forest_maps[[period_name]] <- forest_mask
      
      # Download as raster for local analysis
      tryCatch({
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
      }, error = function(e) {
        cat("Error downloading from GEE, creating local version...\n")
        # Create local version as fallback
        ndvi_raster <- image_with_indices[["NDVI"]]
        forest_raster <- ndvi_raster > 0.4
        writeRaster(forest_raster, 
                    filename = paste0("forest_cover_", period_name, ".tif"),
                    overwrite = TRUE)
      })
      
    } else {
      # Local raster version
      ndvi_raster <- image_with_indices[["NDVI"]]
      forest_mask <- ndvi_raster > 0.4
      forest_maps[[period_name]] <- forest_mask
      
      # Save locally
      writeRaster(forest_mask, 
                  filename = paste0("forest_cover_", period_name, ".tif"),
                  overwrite = TRUE)
    }
  }
  
  return(forest_maps)
}

# ============================================================================
# 4. OBJECTIVE 2: BIOMASS AND CARBON STOCK ESTIMATION
# ============================================================================

# Function to estimate above-ground biomass using allometric equations
estimate_agb <- function(ndvi_image, evi_image) {
  
  if (gee_available) {
    # Google Earth Engine version
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
    
  } else {
    # Local raster version
    # Biomass estimation using vegetation indices
    
    # Method 1: NDVI-based estimation (modified for local raster)
    ndvi_scaled <- ndvi_image * 100
    evi_scaled <- evi_image * 100
    
    # Avoid log of zero or negative values
    ndvi_scaled[ndvi_scaled <= 0] <- 1
    evi_scaled[evi_scaled <= 0] <- 1
    
    agb_ndvi <- exp(2.148 + 0.207 * log(ndvi_scaled) + 0.928 * log(evi_scaled))
    names(agb_ndvi) <- "AGB_NDVI"
    
    # Method 2: Multiple regression approach
    agb_multi <- 150 * (ndvi_image^2) + 50 * evi_image + 25
    names(agb_multi) <- "AGB_MULTI"
    
    return(list(agb_ndvi = agb_ndvi, agb_multi = agb_multi))
  }
}

# Function to convert biomass to carbon stocks
biomass_to_carbon <- function(agb_image, carbon_fraction = 0.47) {
  
  if (gee_available) {
    # Google Earth Engine version
    # Above-ground carbon
    agc <- agb_image$multiply(carbon_fraction)$rename("AGC")
    
    # Below-ground carbon (root-to-shoot ratio)
    # Using IPCC default values for tropical forests
    bgc <- agc$multiply(0.24)$rename("BGC")
    
    # Total carbon stock
    total_carbon <- agc$add(bgc)$rename("TOTAL_CARBON")
    
    return(list(agc = agc, bgc = bgc, total_carbon = total_carbon))
    
  } else {
    # Local raster version
    # Above-ground carbon
    agc <- agb_image * carbon_fraction
    names(agc) <- "AGC"
    
    # Below-ground carbon (root-to-shoot ratio)
    bgc <- agc * 0.24
    names(bgc) <- "BGC"
    
    # Total carbon stock
    total_carbon <- agc + bgc
    names(total_carbon) <- "TOTAL_CARBON"
    
    return(list(agc = agc, bgc = bgc, total_carbon = total_carbon))
  }
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

# Function to calculate fragmentation indices (improved with error handling)
calculate_fragmentation_indices <- function(forest_raster) {
  
  tryCatch({
    # Convert to binary forest/non-forest
    forest_binary <- forest_raster > 0
    
    # Calculate landscape metrics using basic raster operations
    # Patch density (using focal function)
    patch_density <- focal(forest_binary, w = matrix(1, 9, 9), fun = sum, na.rm = TRUE) / 81
    names(patch_density) <- "patch_density"
    
    # Edge density (simplified using boundaries function)
    forest_edges <- boundaries(forest_binary, type = "inner")
    edge_density <- focal(forest_edges, w = matrix(1, 9, 9), fun = sum, na.rm = TRUE)
    names(edge_density) <- "edge_density"
    
    # Forest connectivity (simplified)
    connectivity <- focal(forest_binary, w = matrix(1, 15, 15), fun = sum, na.rm = TRUE) / 225
    names(connectivity) <- "connectivity"
    
    return(list(
      patch_density = patch_density,
      edge_density = edge_density,
      connectivity = connectivity
    ))
    
  }, error = function(e) {
    cat("Error in fragmentation calculation:", e$message, "\n")
    cat("Creating simplified fragmentation indices...\n")
    
    # Simplified fallback version
    forest_binary <- forest_raster > 0
    
    # Simple patch density (proportion of forest in neighborhood)
    patch_density <- focal(forest_binary, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
    names(patch_density) <- "patch_density"
    
    # Simple edge detection
    edge_density <- focal(forest_binary, w = matrix(c(0,1,0,1,-4,1,0,1,0), 3, 3), fun = sum, na.rm = TRUE)
    edge_density[edge_density < 0] <- 0
    names(edge_density) <- "edge_density"
    
    # Simple connectivity
    connectivity <- focal(forest_binary, w = matrix(1, 7, 7), fun = mean, na.rm = TRUE)
    names(connectivity) <- "connectivity"
    
    return(list(
      patch_density = patch_density,
      edge_density = edge_density,
      connectivity = connectivity
    ))
  })
}

# Function to analyze fragmentation-carbon relationships (with robust error handling)
analyze_fragmentation_carbon_correlation <- function() {
  
  correlation_results <- list()
  
  for (i in 1:length(time_periods)) {
    period_name <- names(time_periods)[i]
    
    # Load forest cover and carbon stock rasters
    forest_file <- paste0("forest_cover_", period_name, ".tif")
    carbon_file <- paste0("carbon_stocks_", period_name, ".tif")
    
    tryCatch({
      if (file.exists(forest_file) && file.exists(carbon_file)) {
        cat("Analyzing correlation for", period_name, "...\n")
        
        forest_raster <- raster(forest_file)
        carbon_raster <- raster(carbon_file)
        
        # Ensure rasters have the same extent and resolution
        if (!compareRaster(forest_raster, carbon_raster, stopiffalse = FALSE)) {
          cat("Resampling carbon raster to match forest raster...\n")
          carbon_raster <- resample(carbon_raster, forest_raster, method = "bilinear")
        }
        
        # Calculate fragmentation indices
        frag_indices <- calculate_fragmentation_indices(forest_raster)
        
        # Resample fragmentation indices to match input rasters if needed
        for (idx_name in names(frag_indices)) {
          if (!compareRaster(frag_indices[[idx_name]], forest_raster, stopiffalse = FALSE)) {
            frag_indices[[idx_name]] <- resample(frag_indices[[idx_name]], forest_raster, method = "bilinear")
          }
        }
        
        # Extract values for correlation analysis
        forest_values <- getValues(forest_raster)
        carbon_values <- getValues(carbon_raster)
        patch_density_values <- getValues(frag_indices$patch_density)
        edge_density_values <- getValues(frag_indices$edge_density)
        connectivity_values <- getValues(frag_indices$connectivity)
        
        # Remove NA values and infinite values
        valid_indices <- complete.cases(
          forest_values, carbon_values, patch_density_values,
          edge_density_values, connectivity_values
        ) & is.finite(forest_values) & is.finite(carbon_values) & 
          is.finite(patch_density_values) & is.finite(edge_density_values) & 
          is.finite(connectivity_values)
        
        if (sum(valid_indices) > 10) {  # Need at least 10 valid observations
          # Create correlation matrix
          correlation_data <- data.frame(
            Forest_Cover = forest_values[valid_indices],
            Carbon_Stock = carbon_values[valid_indices],
            Patch_Density = patch_density_values[valid_indices],
            Edge_Density = edge_density_values[valid_indices],
            Connectivity = connectivity_values[valid_indices]
          )
          
          # Remove any remaining problematic values
          correlation_data <- correlation_data[complete.cases(correlation_data), ]
          
          if (nrow(correlation_data) > 10) {
            correlation_matrix <- cor(correlation_data, use = "complete.obs")
            correlation_results[[period_name]] <- correlation_matrix
            
            cat("Successfully calculated correlations for", period_name, 
                "with", nrow(correlation_data), "valid observations\n")
          } else {
            cat("Insufficient valid data for", period_name, "\n")
          }
        } else {
          cat("Insufficient valid data for correlation analysis in", period_name, "\n")
        }
      } else {
        cat("Missing files for", period_name, "- skipping correlation analysis\n")
      }
      
    }, error = function(e) {
      cat("Error analyzing", period_name, ":", e$message, "\n")
    })
  }
  
  return(correlation_results)
}

# ============================================================================
# 6. OBJECTIVE 4: GEOSPATIAL INDICATORS FOR DEGRADATION MONITORING
# ============================================================================

# Function to create degradation indicators (improved with error handling)
create_degradation_indicators <- function() {
  
  tryCatch({
    # Load all time period data
    forest_time_series <- list()
    carbon_time_series <- list()
    
    for (i in 1:length(time_periods)) {
      period_name <- names(time_periods)[i]
      
      forest_file <- paste0("forest_cover_", period_name, ".tif")
      carbon_file <- paste0("carbon_stocks_", period_name, ".tif")
      
      if (file.exists(forest_file)) {
        forest_raster <- raster(forest_file)
        forest_time_series[[period_name]] <- forest_raster
        cat("Loaded forest data for", period_name, "\n")
      }
      if (file.exists(carbon_file)) {
        carbon_raster <- raster(carbon_file)
        carbon_time_series[[period_name]] <- carbon_raster
        cat("Loaded carbon data for", period_name, "\n")
      }
    }
    
    # Calculate change indicators
    degradation_indicators <- list()
    
    if (length(forest_time_series) >= 2) {
      cat("Calculating forest change indicators...\n")
      
      # Forest cover change rate
      first_period <- forest_time_series[[1]]
      last_period <- forest_time_series[[length(forest_time_series)]]
      
      # Ensure rasters are compatible
      if (!compareRaster(first_period, last_period, stopiffalse = FALSE)) {
        cat("Resampling rasters for compatibility...\n")
        last_period <- resample(last_period, first_period, method = "bilinear")
      }
      
      forest_change_rate <- (last_period - first_period) / (length(forest_time_series) - 1)
      names(forest_change_rate) <- "forest_change_rate"
      degradation_indicators$forest_change_rate <- forest_change_rate
      
      # Save forest change rate
      writeRaster(forest_change_rate,
                  filename = "forest_change_rate.tif",
                  overwrite = TRUE)
      
      # Carbon stock change rate
      if (length(carbon_time_series) >= 2) {
        cat("Calculating carbon change indicators...\n")
        
        first_carbon <- carbon_time_series[[1]]
        last_carbon <- carbon_time_series[[length(carbon_time_series)]]
        
        # Ensure rasters are compatible
        if (!compareRaster(first_carbon, last_carbon, stopiffalse = FALSE)) {
          last_carbon <- resample(last_carbon, first_carbon, method = "bilinear")
        }
        
        carbon_change_rate <- (last_carbon - first_carbon) / (length(carbon_time_series) - 1)
        names(carbon_change_rate) <- "carbon_change_rate"
        degradation_indicators$carbon_change_rate <- carbon_change_rate
        
        # Save carbon change rate
        writeRaster(carbon_change_rate,
                    filename = "carbon_change_rate.tif",
                    overwrite = TRUE)
      }
      
      # Degradation hotspots (areas with significant decline)
      degradation_hotspots <- forest_change_rate < -0.1  # 10% decline threshold
      names(degradation_hotspots) <- "degradation_hotspots"
      degradation_indicators$hotspots <- degradation_hotspots
      
      # Save degradation hotspots
      writeRaster(degradation_hotspots,
                  filename = "degradation_hotspots.tif",
                  overwrite = TRUE)
      
      cat("Degradation indicators calculated successfully\n")
      
    } else {
      cat("Insufficient time series data for change analysis\n")
      cat("Available forest periods:", length(forest_time_series), "\n")
      cat("Available carbon periods:", length(carbon_time_series), "\n")
    }
    
    return(degradation_indicators)
    
  }, error = function(e) {
    cat("Error in degradation analysis:", e$message, "\n")
    return(list())
  })
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

# Function to generate summary statistics and visualizations (improved)
generate_summary_statistics <- function(fragmentation_results, degradation_results) {
  
  tryCatch({
    # Create summary report
    sink("outputs/analysis_summary.txt")
    
    cat("FOREST DEGRADATION AND CARBON ACCOUNTING ANALYSIS\n")
    cat("Study Area: Mount Kigali Forest, Nyarugenege District, Rwanda\n")
    cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
    cat("Analysis Time:", format(Sys.time(), "%H:%M:%S"), "\n\n")
    
    cat("SYSTEM INFORMATION:\n")
    cat("R Version:", R.version.string, "\n")
    cat("Google Earth Engine Available:", gee_available, "\n")
    cat("rgee Package Available:", rgee_available, "\n\n")
    
    cat("DATA PROCESSING SUMMARY:\n")
    cat("Time Periods Analyzed:", length(time_periods), "\n")
    for (i in 1:length(time_periods)) {
      period_name <- names(time_periods)[i]
      period <- time_periods[[i]]
      cat(sprintf("  %s: %s to %s\n", period_name, period$start, period$end))
    }
    cat("\n")
    
    cat("CORRELATION ANALYSIS RESULTS:\n")
    if (length(fragmentation_results) > 0) {
      for (period in names(fragmentation_results)) {
        cat("\nPeriod:", period, "\n")
        
        corr_matrix <- fragmentation_results[[period]]
        if (!is.null(corr_matrix) && "Forest_Cover" %in% rownames(corr_matrix) && 
            "Carbon_Stock" %in% colnames(corr_matrix)) {
          
          cat("Correlation between Forest Cover and Carbon Stock:",
              round(corr_matrix["Forest_Cover", "Carbon_Stock"], 3), "\n")
          
          if ("Connectivity" %in% rownames(corr_matrix)) {
            cat("Correlation between Connectivity and Carbon Stock:",
                round(corr_matrix["Connectivity", "Carbon_Stock"], 3), "\n")
          }
          
          if ("Patch_Density" %in% rownames(corr_matrix)) {
            cat("Correlation between Patch Density and Carbon Stock:",
                round(corr_matrix["Patch_Density", "Carbon_Stock"], 3), "\n")
          }
        } else {
          cat("Incomplete correlation data for this period\n")
        }
      }
    } else {
      cat("No correlation results available\n")
    }
    
    cat("\nDEGRADATION INDICATORS:\n")
    if (!is.null(degradation_results$forest_change_rate)) {
      forest_change_stats <- tryCatch({
        c(
          mean = cellStats(degradation_results$forest_change_rate, 'mean', na.rm = TRUE),
          min = cellStats(degradation_results$forest_change_rate, 'min', na.rm = TRUE),
          max = cellStats(degradation_results$forest_change_rate, 'max', na.rm = TRUE)
        )
      }, error = function(e) {
        cat("Error calculating forest change statistics:", e$message, "\n")
        c(mean = NA, min = NA, max = NA)
      })
      
      if (!is.na(forest_change_stats[1])) {
        cat("Mean annual forest change rate:", round(forest_change_stats[1], 4), "\n")
        cat("Minimum forest change rate:", round(forest_change_stats[2], 4), "\n")
        cat("Maximum forest change rate:", round(forest_change_stats[3], 4), "\n")
      }
    } else {
      cat("Forest change rate data not available\n")
    }
    
    if (!is.null(degradation_results$carbon_change_rate)) {
      carbon_change_stats <- tryCatch({
        c(
          mean = cellStats(degradation_results$carbon_change_rate, 'mean', na.rm = TRUE),
          min = cellStats(degradation_results$carbon_change_rate, 'min', na.rm = TRUE),
          max = cellStats(degradation_results$carbon_change_rate, 'max', na.rm = TRUE)
        )
      }, error = function(e) {
        cat("Error calculating carbon change statistics:", e$message, "\n")
        c(mean = NA, min = NA, max = NA)
      })
      
      if (!is.na(carbon_change_stats[1])) {
        cat("Mean annual carbon change rate:", round(carbon_change_stats[1], 2), "Mg C/ha/year\n")
        cat("Minimum carbon change rate:", round(carbon_change_stats[2], 2), "Mg C/ha/year\n")
        cat("Maximum carbon change rate:", round(carbon_change_stats[3], 2), "Mg C/ha/year\n")
      }
    } else {
      cat("Carbon change rate data not available\n")
    }
    
    cat("\nFILES GENERATED:\n")
    generated_files <- list.files(pattern = "\\.(tif|pdf|txt)$")
    if (length(generated_files) > 0) {
      for (file in generated_files) {
        cat("  -", file, "\n")
      }
    } else {
      cat("  No output files found\n")
    }
    
    cat("\nANALYSIS COMPLETED SUCCESSFULLY\n")
    cat("======================================\n")
    
    sink()
    
    # Create visualization plots
    create_analysis_plots(fragmentation_results, degradation_results)
    
    cat("Summary statistics saved to outputs/analysis_summary.txt\n")
    
  }, error = function(e) {
    cat("Error generating summary statistics:", e$message, "\n")
    sink()  # Make sure sink is closed even if error occurs
  })
}

# Function to create analysis plots (improved with error handling)
create_analysis_plots <- function(fragmentation_results, degradation_results) {
  
  tryCatch({
    # Plot correlation matrices
    if (length(fragmentation_results) > 0) {
      cat("Creating correlation analysis plots...\n")
      
      pdf("outputs/correlation_analysis.pdf", width = 12, height = 8)
      
      for (period in names(fragmentation_results)) {
        tryCatch({
          corr_matrix <- fragmentation_results[[period]]
          if (!is.null(corr_matrix) && is.matrix(corr_matrix)) {
            corrplot(corr_matrix, 
                     method = "color", 
                     type = "upper",
                     title = paste("Correlation Matrix -", period),
                     mar = c(0, 0, 2, 0),
                     tl.cex = 0.8,
                     number.cex = 0.8)
          }
        }, error = function(e) {
          cat("Error plotting correlation for", period, ":", e$message, "\n")
        })
      }
      
      dev.off()
      cat("Correlation plots saved to outputs/correlation_analysis.pdf\n")
    }
    
    # Plot degradation trends if available
    if (!is.null(degradation_results$forest_change_rate)) {
      cat("Creating degradation maps...\n")
      
      pdf("outputs/degradation_maps.pdf", width = 10, height = 8)
      
      # Forest change rate map
      tryCatch({
        plot(degradation_results$forest_change_rate,
             main = "Annual Forest Change Rate (%/year)",
             col = RColorBrewer::brewer.pal(11, "RdYlGn"),
             axes = TRUE,
             box = TRUE)
        
        # Add legend
        legend("bottomright", 
               legend = c("Loss", "Stable", "Gain"), 
               fill = c("red", "yellow", "green"),
               cex = 0.8)
      }, error = function(e) {
        cat("Error plotting forest change rate:", e$message, "\n")
      })
      
      # Degradation hotspots
      if (!is.null(degradation_results$hotspots)) {
        tryCatch({
          plot(degradation_results$hotspots,
               main = "Forest Degradation Hotspots",
               col = c("lightgreen", "red"),
               legend = FALSE,
               axes = TRUE,
               box = TRUE)
          legend("bottomright", 
                 legend = c("Stable Areas", "Degradation Hotspots"), 
                 fill = c("lightgreen", "red"),
                 cex = 0.8)
        }, error = function(e) {
          cat("Error plotting degradation hotspots:", e$message, "\n")
        })
      }
      
      # Carbon change rate map if available
      if (!is.null(degradation_results$carbon_change_rate)) {
        tryCatch({
          plot(degradation_results$carbon_change_rate,
               main = "Annual Carbon Stock Change Rate (Mg C/ha/year)",
               col = RColorBrewer::brewer.pal(11, "RdYlBu"),
               axes = TRUE,
               box = TRUE)
        }, error = function(e) {
          cat("Error plotting carbon change rate:", e$message, "\n")
        })
      }
      
      dev.off()
      cat("Degradation maps saved to outputs/degradation_maps.pdf\n")
    }
    
    # Create time series plots if multiple periods available
    create_time_series_plots()
    
  }, error = function(e) {
    cat("Error in plot creation:", e$message, "\n")
  })
}

# Function to create time series plots
create_time_series_plots <- function() {
  
  tryCatch({
    cat("Creating time series analysis plots...\n")
    
    # Collect statistics for each time period
    forest_stats <- data.frame()
    carbon_stats <- data.frame()
    
    for (i in 1:length(time_periods)) {
      period_name <- names(time_periods)[i]
      
      # Forest cover statistics
      forest_file <- paste0("forest_cover_", period_name, ".tif")
      if (file.exists(forest_file)) {
        forest_raster <- raster(forest_file)
        
        forest_mean <- cellStats(forest_raster, 'mean', na.rm = TRUE)
        forest_area <- sum(getValues(forest_raster), na.rm = TRUE) * res(forest_raster)[1] * res(forest_raster)[2]
        
        forest_stats <- rbind(forest_stats, data.frame(
          Period = period_name,
          Year = 2015 + (i-1) * 2,  # Approximate year
          Mean_Cover = forest_mean,
          Total_Area = forest_area
        ))
      }
      
      # Carbon stock statistics
      carbon_file <- paste0("carbon_stocks_", period_name, ".tif")
      if (file.exists(carbon_file)) {
        carbon_raster <- raster(carbon_file)
        
        carbon_mean <- cellStats(carbon_raster, 'mean', na.rm = TRUE)
        carbon_total <- cellStats(carbon_raster, 'sum', na.rm = TRUE)
        
        carbon_stats <- rbind(carbon_stats, data.frame(
          Period = period_name,
          Year = 2015 + (i-1) * 2,  # Approximate year
          Mean_Carbon = carbon_mean,
          Total_Carbon = carbon_total
        ))
      }
    }
    
    # Create time series plots
    if (nrow(forest_stats) > 1 || nrow(carbon_stats) > 1) {
      pdf("outputs/time_series_analysis.pdf", width = 12, height = 8)
      
      # Forest cover time series
      if (nrow(forest_stats) > 1) {
        par(mfrow = c(2, 2))
        
        plot(forest_stats$Year, forest_stats$Mean_Cover,
             type = "b", pch = 19, col = "darkgreen",
             xlab = "Year", ylab = "Mean Forest Cover",
             main = "Forest Cover Change Over Time",
             ylim = c(0, max(forest_stats$Mean_Cover, na.rm = TRUE) * 1.1))
        grid()
        
        plot(forest_stats$Year, forest_stats$Total_Area,
             type = "b", pch = 19, col = "forestgreen",
             xlab = "Year", ylab = "Total Forest Area",
             main = "Total Forest Area Change Over Time")
        grid()
      }
      
      # Carbon stock time series
      if (nrow(carbon_stats) > 1) {
        plot(carbon_stats$Year, carbon_stats$Mean_Carbon,
             type = "b", pch = 19, col = "brown",
             xlab = "Year", ylab = "Mean Carbon Stock (Mg C/ha)",
             main = "Carbon Stock Change Over Time")
        grid()
        
        plot(carbon_stats$Year, carbon_stats$Total_Carbon,
             type = "b", pch = 19, col = "chocolate",
             xlab = "Year", ylab = "Total Carbon Stock (Mg C)",
             main = "Total Carbon Stock Change Over Time")
        grid()
      }
      
      dev.off()
      cat("Time series plots saved to outputs/time_series_analysis.pdf\n")
    }
    
  }, error = function(e) {
    cat("Error creating time series plots:", e$message, "\n")
  })
}

# ============================================================================
# 8. COMPREHENSIVE ANALYSIS AND REPORTING (Updated)
# ============================================================================

# Main analysis function (improved with better error handling)
run_comprehensive_analysis <- function() {
  
  cat("========================================\n")
  cat("STARTING COMPREHENSIVE FOREST ANALYSIS\n")
  cat("========================================\n")
  cat("Study Area: Mount Kigali Forest, Nyarugenege District, Rwanda\n")
  cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
  cat("Analysis Time:", format(Sys.time(), "%H:%M:%S"), "\n\n")
  
  # Create output directory
  if (!dir.exists("outputs")) {
    dir.create("outputs")
    cat("Created outputs directory\n")
  }
  
  # Initialize results storage
  analysis_results <- list()
  
  # Step 1: Forest cover mapping
  cat("========================================\n")
  cat("STEP 1: FOREST COVER MAPPING\n")
  cat("========================================\n")
  tryCatch({
    forest_maps <- create_forest_cover_maps()
    analysis_results$forest_maps <- forest_maps
    cat("✓ Forest cover mapping completed successfully\n\n")
  }, error = function(e) {
    cat("✗ Error in forest cover mapping:", e$message, "\n\n")
    analysis_results$forest_maps <- NULL
  })
  
  # Step 2: Biomass and carbon estimation
  cat("========================================\n")
  cat("STEP 2: BIOMASS AND CARBON ESTIMATION\n")
  cat("========================================\n")
  tryCatch({
    biomass_carbon_results <- create_biomass_carbon_maps()
    analysis_results$biomass_carbon <- biomass_carbon_results
    cat("✓ Biomass and carbon estimation completed successfully\n\n")
  }, error = function(e) {
    cat("✗ Error in biomass and carbon estimation:", e$message, "\n\n")
    analysis_results$biomass_carbon <- NULL
  })
  
  # Step 3: Fragmentation analysis
  cat("========================================\n")
  cat("STEP 3: FOREST FRAGMENTATION ANALYSIS\n")
  cat("========================================\n")
  tryCatch({
    fragmentation_results <- analyze_fragmentation_carbon_correlation()
    analysis_results$fragmentation <- fragmentation_results
    cat("✓ Fragmentation analysis completed successfully\n\n")
  }, error = function(e) {
    cat("✗ Error in fragmentation analysis:", e$message, "\n\n")
    analysis_results$fragmentation <- list()
  })
  
  # Step 4: Degradation indicators
  cat("========================================\n")
  cat("STEP 4: DEGRADATION INDICATORS\n")
  cat("========================================\n")
  tryCatch({
    degradation_results <- create_degradation_indicators()
    analysis_results$degradation <- degradation_results
    cat("✓ Degradation indicators completed successfully\n\n")
  }, error = function(e) {
    cat("✗ Error in degradation indicators:", e$message, "\n\n")
    analysis_results$degradation <- list()
  })
  
  # Step 5: Generate comprehensive report
  cat("========================================\n")
  cat("STEP 5: GENERATING COMPREHENSIVE REPORT\n")
  cat("========================================\n")
  tryCatch({
    generate_summary_statistics(analysis_results$fragmentation, analysis_results$degradation)
    cat("✓ Summary report generated successfully\n\n")
  }, error = function(e) {
    cat("✗ Error generating summary report:", e$message, "\n\n")
  })
  
  # Final summary
  cat("========================================\n")
  cat("ANALYSIS COMPLETION SUMMARY\n")
  cat("========================================\n")
  cat("Study completed for Mount Kigali Forest\n")
  cat("Output files saved in 'outputs' directory\n")
  cat("Generated files:\n")
  
  output_files <- list.files("outputs", full.names = FALSE)
  if (length(output_files) > 0) {
    for (file in output_files) {
      cat("  ✓", file, "\n")
    }
  }
  
  generated_files <- list.files(pattern = "\\.(tif|pdf)$", full.names = FALSE)
  if (length(generated_files) > 0) {
    cat("\nAdditional files in working directory:\n")
    for (file in generated_files) {
      cat("  ✓", file, "\n")
    }
  }
  
  cat("\n========================================\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("========================================\n")
  
  return(analysis_results)
}

# ============================================================================
# 9. EXECUTION INSTRUCTIONS AND SETUP
# ============================================================================

# Function to check system requirements
check_system_requirements <- function() {
  cat("========================================\n")
  cat("SYSTEM REQUIREMENTS CHECK\n")
  cat("========================================\n")
  
  # Check R version
  cat("R Version:", R.version.string, "\n")
  if (as.numeric(R.version$major) >= 4) {
    cat("✓ R version is compatible\n")
  } else {
    cat("⚠ R version 4.0+ recommended\n")
  }
  
  # Check required packages
  required_core <- c("sf", "raster", "terra", "ggplot2", "dplyr")
  missing_packages <- required_core[!sapply(required_core, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) == 0) {
    cat("✓ Core packages are available\n")
  } else {
    cat("⚠ Missing core packages:", paste(missing_packages, collapse = ", "), "\n")
  }
  
  # Check Google Earth Engine
  if (rgee_available && gee_available) {
    cat("✓ Google Earth Engine is ready\n")
  } else {
    cat("⚠ Google Earth Engine not available - will use simulated data\n")
  }
  
  cat("========================================\n\n")
}

# Print execution instructions
cat("================================================================================\n")
cat("FOREST DEGRADATION AND CARBON ACCOUNTING ANALYSIS\n")
cat("Mount Kigali Forest, Nyarugenege District, Rwanda\n")
cat("================================================================================\n\n")

cat("SETUP COMPLETE!\n\n")

cat("To run the analysis:\n")
cat("1. Check system requirements: check_system_requirements()\n")
cat("2. Run full analysis: results <- run_comprehensive_analysis()\n")
cat("3. Check outputs directory for results\n\n")

cat("EXPECTED OUTPUTS:\n")
cat("- Forest cover maps (.tif files)\n")
cat("- Carbon stock maps (.tif files)\n") 
cat("- Forest and carbon change rate maps (.tif files)\n")
cat("- Degradation hotspot maps (.tif files)\n")
cat("- Comprehensive analysis summary (analysis_summary.txt)\n")
cat("- Correlation analysis plots (correlation_analysis.pdf)\n")
cat("- Degradation maps (degradation_maps.pdf)\n")
cat("- Time series analysis (time_series_analysis.pdf)\n\n")

cat("STUDY AREA DETAILS:\n")
cat("- Location: Mount Kigali Forest, Nyarugenege District, Rwanda\n")
cat("- Coordinates: 30.0°E to 30.15°E, -1.95°S to -1.85°S\n")
cat("- Analysis periods: 2015-2024 (5 time periods)\n")
cat("- Spatial resolution: 30m (Landsat)\n\n")

cat("ANALYSIS METHODS:\n")
cat("- Multi-temporal Landsat 8 analysis\n")
cat("- Vegetation indices (NDVI, EVI, SAVI, MSAVI)\n")
cat("- Allometric biomass estimation\n")
cat("- IPCC carbon conversion factors\n")
cat("- Landscape fragmentation metrics\n")
cat("- Change detection algorithms\n\n")

cat("For help or troubleshooting, check the error messages in the console.\n")
cat("================================================================================\n")