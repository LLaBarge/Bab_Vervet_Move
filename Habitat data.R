# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(mgcv)
library(terra)
library(sf)
library(ggplot2)
library(viridis)

# Set working directory
setwd("C:/Users/lrlab/PPP scan points")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

# Load vegetation plot data
veg_plots <- read_excel("vegplot_phenology.xlsx", sheet = 1)

# Load phenology data
phenology <- read_excel("vegplot_phenology.xlsx", sheet = "phenology")

# IMPORTANT: Replace all NA values with 0 at the beginning
veg_plots[is.na(veg_plots)] <- 0
phenology[is.na(phenology)] <- 0

cat("Vegetation plots loaded:", nrow(veg_plots), "\n")
cat("Phenology observations loaded:", nrow(phenology), "\n")

# ============================================================================
# 2. CLEAN AND FILTER PHENOLOGY DATA
# ============================================================================

# Standardize month names (convert to title case)
phenology$Month <- tools::toTitleCase(tolower(trimws(phenology$Month)))

# Convert month names to numbers
month_mapping <- c(
  "January" = 1, "February" = 2, "March" = 3, "April" = 4,
  "May" = 5, "June" = 6, "July" = 7, "August" = 8,
  "September" = 9, "October" = 10, "November" = 11, "December" = 12
)
phenology$Month_Num <- month_mapping[phenology$Month]

# Filter out years 2012, 2013, and 2019+
phenology_filtered <- phenology %>%
  filter(Year >= 2014 & Year <= 2018)

cat("Years kept:", unique(phenology_filtered$Year), "\n")
cat("Phenology observations after filtering:", nrow(phenology_filtered), "\n")

# Convert fruit and seed columns to numeric (handle any text entries)
phenology_filtered$`No. fruit on tree` <- as.numeric(phenology_filtered$`No. fruit on tree`)
phenology_filtered$`No. seed pods on tree` <- as.numeric(phenology_filtered$`No. seed pods on tree`)
phenology_filtered[is.na(phenology_filtered)] <- 0

# ============================================================================
# 3. IDENTIFY SPECIES COLUMNS AND PREPARE VEG DATA
# ============================================================================

# Identify species columns
species_cols <- names(veg_plots)[grepl("^[A-Z]\\.", names(veg_plots))]
cat("\nSpecies columns found:", length(species_cols), "\n")

# Calculate total tree count and species richness per plot
veg_plots$total_trees <- rowSums(veg_plots[, species_cols])
veg_plots$species_richness <- rowSums(veg_plots[, species_cols] > 0)

# ============================================================================
# 4. CONVERT TO SPATIAL AND GET UTM COORDINATES
# ============================================================================

cat("\n=== Converting to Spatial Format ===\n")

# Convert to sf object (WGS84)
veg_sf <- st_as_sf(veg_plots, coords = c("Long", "Lat"), crs = 4326)

# Convert to UTM Zone 36S (adjust if different location)
utm_crs <- 32736  
veg_utm <- st_transform(veg_sf, crs = utm_crs)

# Extract UTM coordinates
coords <- st_coordinates(veg_utm)
veg_plots$UTM_X <- coords[, 1]
veg_plots$UTM_Y <- coords[, 2]

# ============================================================================
# 5. CREATE PREDICTION GRID
# ============================================================================

cat("\n=== Creating Prediction Grid ===\n")


x_range <- range(veg_plots$UTM_X)
y_range <- range(veg_plots$UTM_Y)

x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05

resolution <- 10  # meters

prediction_grid <- expand.grid(
  UTM_X = seq(x_range[1] - x_buffer, x_range[2] + x_buffer, by = resolution),
  UTM_Y = seq(y_range[1] - y_buffer, y_range[2] + y_buffer, by = resolution)
)

cat("Prediction grid created:", nrow(prediction_grid), "cells\n")



# ============================================================================
# 6. BUILD STATIC GAM MODELS (Canopy metrics - don't change much)
# ============================================================================
#rename the canopy cover variable
veg_plots$CC<-veg_plots$`Canopy.Cover(%)`

# Model 1: Canopy Cover
cat("Fitting Canopy Cover model...\n")
gam_canopy_cover <- gam(CC ~ s(UTM_X, UTM_Y, k = 30) + 
                          s(species_richness, k = 5),
                        data = veg_plots,
                        method = "REML"
)
CH<-veg_plots$`canopy.height(m)`

# Model 2: Canopy Height
cat("Fitting Canopy Height model...\n")
gam_canopy_height <- gam(CH ~ s(UTM_X, UTM_Y, k = 30) + 
                           s(CC, k = 5),
                         data = veg_plots,
                         method = "REML"
)

HV<-veg_plots$`Horizontal.vis(m)`
# Model 3: Horizontal Visibility
cat("Fitting Horizontal Visibility model...\n")
gam_horizontal_vis <- gam(HV ~ s(UTM_X, UTM_Y, k = 30) + 
                            s(CC, k = 5) +
                            s(species_richness, k = 5),
                          data = veg_plots,
                          method = "REML"
)

# ============================================================================
# 7. PREDICT STATIC VARIABLES
# ============================================================================

cat("\n=== Predicting Static Variables ===\n")

# Add mean species richness for predictions
prediction_grid$species_richness <- mean(veg_plots$species_richness)

# Predict canopy cover
prediction_grid$canopy_cover_pred <- predict(gam_canopy_cover, 
                                             newdata = prediction_grid, 
                                             type = "response")
prediction_grid$`Canopy.Cover(%)` <- prediction_grid$canopy_cover_pred

# Predict canopy height
prediction_grid$canopy_height_pred <- predict(gam_canopy_height, 
                                              newdata = prediction_grid, 
                                              type = "response")
prediction_grid$`canopy.height(m)` <- prediction_grid$canopy_height_pred

# Predict horizontal visibility
prediction_grid$horizontal_vis_pred <- predict(gam_horizontal_vis, 
                                               newdata = prediction_grid, 
                                               type = "response")

# ============================================================================
# 8. CREATE STATIC RASTERS
# ============================================================================

cat("\n=== Creating Static Rasters ===\n")

# Create raster template
raster_template <- rast(
  xmin = min(prediction_grid$UTM_X), 
  xmax = max(prediction_grid$UTM_X),
  ymin = min(prediction_grid$UTM_Y), 
  ymax = max(prediction_grid$UTM_Y),
  resolution = resolution,
  crs = paste0("EPSG:", utm_crs)
)

# Convert to spatial vector for rasterization
pred_vect <- vect(prediction_grid, geom = c("UTM_X", "UTM_Y"), 
                  crs = paste0("EPSG:", utm_crs))

# Create rasters
canopy_cover_raster <- rasterize(pred_vect, raster_template, 
                                 field = "Canopy.Cover(%)")
canopy_height_raster <- rasterize(pred_vect, raster_template, 
                                  field = "canopy.height(m)")
horizontal_vis_raster <- rasterize(pred_vect, raster_template, 
                                   field = "horizontal_vis_pred")

# Save static rasters
dir.create("raster_outputs", showWarnings = FALSE)

writeRaster(canopy_cover_raster, "raster_outputs/canopy_cover.tif", overwrite = TRUE)
writeRaster(canopy_height_raster, "raster_outputs/canopy_height.tif", overwrite = TRUE)
writeRaster(horizontal_vis_raster, "raster_outputs/horizontal_visibility.tif", overwrite = TRUE)

cat("Static rasters saved\n")

# ============================================================================
# 9. CREATE YEARLY-MONTHLY FRUIT AND SEED RASTERS
# ============================================================================

cat("\n=== Creating Yearly-Monthly Fruit and Seed Rasters ===\n")

# Get unique years and months
years <- sort(unique(phenology_filtered$Year))

# Loop through each year
for (year in years) {
  
  cat("\n### Processing Year:", year, "###\n")
  
  # Filter phenology for this year
  pheno_year <- phenology_filtered %>% filter(Year == year)
  
  # Get available months for this year
  months_in_year <- sort(unique(pheno_year$Month_Num))
  
  # Loop through each month
  for (month_num in months_in_year) {
    
    month_name <- names(month_mapping)[month_mapping == month_num]
    cat("Processing", month_name, year, "...\n")
    
    # Filter phenology for this year-month
    pheno_month <- pheno_year %>% filter(Month_Num == month_num)
    
    # Calculate species productivity for this month
    species_productivity <- pheno_month %>%
      group_by(Species) %>%
      summarise(
        mean_fruit_per_tree = mean(`No. fruit on tree`),
        mean_seed_per_tree = mean(`No. seed pods on tree`),
        n_trees_sampled = n(),
        .groups = "drop"
      )
    
    # Calculate estimated fruit and seed availability per plot
    veg_month <- veg_plots
    veg_month$estimated_fruit <- 0
    veg_month$estimated_seeds <- 0
    
    for (species in species_cols) {
      if (species %in% species_productivity$Species) {
        species_data <- species_productivity %>% filter(Species == species)
        
        veg_month$estimated_fruit <- veg_month$estimated_fruit + 
          (veg_month[[species]] * species_data$mean_fruit_per_tree)
        
        veg_month$estimated_seeds <- veg_month$estimated_seeds + 
          (veg_month[[species]] * species_data$mean_seed_per_tree)
      }
    }
    
    # Only model if we have enough data points with fruit/seeds
    plots_with_fruit <- sum(veg_month$estimated_fruit > 0)
    plots_with_seeds <- sum(veg_month$estimated_seeds > 0)
    
    cat("  Plots with fruit:", plots_with_fruit, 
        "| Plots with seeds:", plots_with_seeds, "\n")
    
    # Add covariates
    veg_month$total_trees <- rowSums(veg_month[, species_cols])
    
    # ========================================================================
    # FRUIT MODEL - More localized (smaller k for spatial smooth)
    # ========================================================================
    
    if (plots_with_fruit >= 20) {
      
      veg_fruit <- veg_month %>% filter(estimated_fruit > 0)
      
      # Use SMALLER k for spatial smooth to capture localized fruiting patterns
      gam_fruit <- gam(
        log(estimated_fruit + 1) ~ s(UTM_X, UTM_Y, k = 15) +  # REDUCED from 30
          s(`Canopy.Cover(%)`, k = 5) +
          s(total_trees, k = 5),
        data = veg_fruit,
        method = "REML"
      )
      
      # Prepare prediction grid with covariates
      pred_grid_month <- prediction_grid
      pred_grid_month$total_trees <- mean(veg_month$total_trees)
      
      # Predict (log scale)
      pred_grid_month$fruit_pred_log <- predict(gam_fruit, 
                                                newdata = pred_grid_month, 
                                                type = "response")
      
      # Back-transform
      pred_grid_month$fruit_availability <- exp(pred_grid_month$fruit_pred_log) - 1
      
      # Create raster
      pred_vect_month <- vect(pred_grid_month, geom = c("UTM_X", "UTM_Y"), 
                              crs = paste0("EPSG:", utm_crs))
      
      fruit_raster_month <- rasterize(pred_vect_month, raster_template, 
                                      field = "fruit_availability")
      
      # Save
      dir.create(paste0("raster_outputs/", year), showWarnings = FALSE)
      filename <- paste0("raster_outputs/", year, "/fruit_", 
                         tolower(month_name), ".tif")
      writeRaster(fruit_raster_month, filename, overwrite = TRUE)
      
      cat("  Saved:", filename, "\n")
    } else {
      cat("  Skipping fruit raster (insufficient data)\n")
    }
    
    # ========================================================================
    # SEED MODEL - More localized (smaller k for spatial smooth)
    # ========================================================================
    
    if (plots_with_seeds >= 20) {
      
      veg_seed <- veg_month %>% filter(estimated_seeds > 0)
      
      # Use SMALLER k for spatial smooth to capture localized seeding patterns
      gam_seed <- gam(
        log(estimated_seeds + 1) ~ s(UTM_X, UTM_Y, k = 15) +  # REDUCED from 30
          s(`Canopy.Cover(%)`, k = 5) +
          s(total_trees, k = 5),
        data = veg_seed,
        method = "REML"
      )
      
      # Prepare prediction grid
      pred_grid_month <- prediction_grid
      pred_grid_month$total_trees <- mean(veg_month$total_trees)
      
      # Predict (log scale)
      pred_grid_month$seed_pred_log <- predict(gam_seed, 
                                               newdata = pred_grid_month, 
                                               type = "response")
      
      # Back-transform
      pred_grid_month$seed_availability <- exp(pred_grid_month$seed_pred_log) - 1
      
      # Create raster
      pred_vect_month <- vect(pred_grid_month, geom = c("UTM_X", "UTM_Y"), 
                              crs = paste0("EPSG:", utm_crs))
      
      seed_raster_month <- rasterize(pred_vect_month, raster_template, 
                                     field = "seed_availability")
      
      # Save
      filename <- paste0("raster_outputs/", year, "/seed_", 
                         tolower(month_name), ".tif")
      writeRaster(seed_raster_month, filename, overwrite = TRUE)
      
      cat("  Saved:", filename, "\n")
    } else {
      cat("  Skipping seed raster (insufficient data)\n")
    }
  }
}

# ============================================================================
# 10. CREATE VISUALIZATIONS
# ============================================================================

cat("\n=== Creating Visualizations ===\n")

# Plot static rasters
png("raster_outputs/canopy_cover_map.png", width = 800, height = 800)
plot(canopy_cover_raster, main = "Canopy Cover (%)", col = viridis(100))
points(veg_plots$UTM_X, veg_plots$UTM_Y, pch = 20, cex = 0.3, col = "red")
dev.off()

png("raster_outputs/canopy_height_map.png", width = 800, height = 800)
plot(canopy_height_raster, main = "Canopy Height (m)", col = viridis(100))
points(veg_plots$UTM_X, veg_plots$UTM_Y, pch = 20, cex = 0.3, col = "red")
dev.off()

png("raster_outputs/horizontal_visibility_map.png", width = 800, height = 800)
plot(horizontal_vis_raster, main = "Horizontal Visibility (m)", col = viridis(100))
points(veg_plots$UTM_X, veg_plots$UTM_Y, pch = 20, cex = 0.3, col = "red")
dev.off()

# Create example visualizations for one year-month
example_files <- list.files("raster_outputs", pattern = "fruit.*\\.tif$", 
                            recursive = TRUE, full.names = TRUE)
if (length(example_files) > 0) {
  example_raster <- rast(example_files[1])
  png("raster_outputs/example_fruit_availability.png", width = 800, height = 800)
  plot(example_raster, main = paste("Example:", basename(example_files[1])), 
       col = viridis(100))
  points(veg_plots$UTM_X, veg_plots$UTM_Y, pch = 20, cex = 0.3, col = "red")
  dev.off()
}
