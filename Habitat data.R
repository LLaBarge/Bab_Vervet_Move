# This version processes predictions in chunks to avoid memory issues

library(readxl)
library(dplyr)
library(mgcv)
library(terra)
library(sf)

setwd("C:/Users/lrlab/PPP scan points")


# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading data...\n")
veg_plots <- read_excel("vegplot_phenology.xlsx", sheet = 1)
phenology <- read_excel("vegplot_phenology.xlsx", sheet = "phenology")

# Replace NA with 0
veg_plots[is.na(veg_plots)] <- 0
phenology[is.na(phenology)] <- 0

# ============================================================================
# 2. FILTER PHENOLOGY
# ============================================================================

phenology$Month <- tools::toTitleCase(tolower(trimws(phenology$Month)))
month_mapping <- c(
  "January" = 1, "February" = 2, "March" = 3, "April" = 4,
  "May" = 5, "June" = 6, "July" = 7, "August" = 8,
  "September" = 9, "October" = 10, "November" = 11, "December" = 12
)
phenology$Month_Num <- month_mapping[phenology$Month]
phenology_filtered <- phenology %>% filter(Year >= 2014 & Year <= 2018)
phenology_filtered$`No. fruit on tree` <- as.numeric(phenology_filtered$`No. fruit on tree`)
phenology_filtered$`No. seed pods on tree` <- as.numeric(phenology_filtered$`No. seed pods on tree`)
phenology_filtered[is.na(phenology_filtered)] <- 0

# ============================================================================
# 3. PREPARE VEG DATA
# ============================================================================

species_cols <- names(veg_plots)[grepl("^[A-Z]\\.", names(veg_plots))]
veg_plots$total_trees <- rowSums(veg_plots[, species_cols])
veg_plots$species_richness <- rowSums(veg_plots[, species_cols] > 0)

# Get UTM coordinates
veg_sf <- st_as_sf(veg_plots, coords = c("Long", "Lat"), crs = 4326)
veg_utm <- st_transform(veg_sf, crs = 32736)
coords <- st_coordinates(veg_utm)
veg_plots$UTM_X <- coords[, 1]
veg_plots$UTM_Y <- coords[, 2]

# ============================================================================
# 4. CREATE RASTER TEMPLATE (NO DATA FRAME GRID)
# ============================================================================

# create the template raster

x_range <- range(veg_plots$UTM_X)
y_range <- range(veg_plots$UTM_Y)
x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05

# 25m resolution
resolution <- 25  # 50 meters - maybe adjust as needed

raster_template <- rast(
  xmin = x_range[1] - x_buffer, 
  xmax = x_range[2] + x_buffer,
  ymin = y_range[1] - y_buffer, 
  ymax = y_range[2] + y_buffer,
  resolution = resolution,
  crs = "EPSG:32736"
)

cat("Raster size:", nrow(raster_template), "x", ncol(raster_template), "\n")
cat("Total cells:", ncell(raster_template), "\n\n")

# ============================================================================
# 5. BUILD STATIC MODELS
# ============================================================================

cat("Building static GAM models...\n")

#rename these variables
veg_plots$CH<-veg_plots$`canopy.height(m)`
veg_plots$CC<-veg_plots$`Canopy.Cover(%)`
veg_plots$HV<-veg_plots$`Horizontal.vis(m)`

gam_canopy_cover <- gam(CC ~ s(UTM_X, UTM_Y, k = 25) + s(species_richness, k = 5),
  data = veg_plots, method = "REML"
)

gam_canopy_height <- gam(CH ~ s(UTM_X, UTM_Y, k = 25) + s(CC, k = 5),
  data = veg_plots, method = "REML"
)

gam_horizontal_vis <- gam(HV ~ s(UTM_X, UTM_Y, k = 25) + 
    s(CC, k = 5) + 
    s(species_richness, k = 5),
  data = veg_plots, method = "REML"
)

# ============================================================================
# 6. PREDICT DIRECTLY TO RASTERS (NO BIG DATA FRAME)
# ============================================================================

dir.create("raster_outputs", showWarnings = FALSE)

# Function to predict to raster in chunks
predict_to_raster <- function(model, template, extra_vars = NULL) {
  
  # Get coordinates from raster
  coords_matrix <- crds(template, df = FALSE)
  
  # Process in chunks to avoid memory issues
  chunk_size <- 50000  # Adjust if needed
  n_cells <- nrow(coords_matrix)
  n_chunks <- ceiling(n_cells / chunk_size)
  
  predictions <- numeric(n_cells)
  
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_cells)
    
    # Create prediction data for this chunk
    chunk_data <- data.frame(
      UTM_X = coords_matrix[start_idx:end_idx, 1],
      UTM_Y = coords_matrix[start_idx:end_idx, 2]
    )
    
    # Add extra variables if provided
    if (!is.null(extra_vars)) {
      for (var_name in names(extra_vars)) {
        chunk_data[[var_name]] <- extra_vars[[var_name]]
      }
    }
    
    # Predict for this chunk
    predictions[start_idx:end_idx] <- predict(model, newdata = chunk_data, 
                                              type = "response")
    
    if (i %% 10 == 0) {
      cat("  Processed chunk", i, "of", n_chunks, "\n")
    }
  }
  
  # Assign values to raster
  values(template) <- predictions
  return(template)
}

# Create static rasters
cat("Canopy cover...\n")
canopy_cover_raster <- predict_to_raster(
  gam_canopy_cover, 
  raster_template, 
  extra_vars = list(species_richness = mean(veg_plots$species_richness))
)
writeRaster(canopy_cover_raster, "raster_outputs/canopy_cover.tif", 
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("Canopy height...\n")
canopy_height_raster <- predict_to_raster(
  gam_canopy_height, 
  raster_template,
  extra_vars = list(`Canopy.Cover(%)` = mean(veg_plots$`Canopy.Cover(%)`))
)
writeRaster(canopy_height_raster, "raster_outputs/canopy_height.tif", 
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("Horizontal visibility...\n")
horizontal_vis_raster <- predict_to_raster(
  gam_horizontal_vis, 
  raster_template,
  extra_vars = list(
    `Canopy.Cover(%)` = mean(veg_plots$`Canopy.Cover(%)`),
    species_richness = mean(veg_plots$species_richness)
  )
)
writeRaster(horizontal_vis_raster, "raster_outputs/horizontal_visibility.tif", 
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# Get canopy cover values for fruit/seed models
canopy_values <- values(canopy_cover_raster)[, 1]

# Clean up
rm(gam_canopy_cover, gam_canopy_height, gam_horizontal_vis)
gc()


# ============================================================================
# 7. YEARLY-MONTHLY FRUIT AND SEED RASTERS
# ============================================================================

cat("=== Creating Monthly Rasters ===\n")

years <- sort(unique(phenology_filtered$Year))

for (year in years) {
  
  cat("\n### Year:", year, "###\n")
  
  pheno_year <- phenology_filtered %>% filter(Year == year)
  months_in_year <- sort(unique(pheno_year$Month_Num))
  
  for (month_num in months_in_year) {
    
    month_name <- names(month_mapping)[month_mapping == month_num]
    cat(month_name, "...\n")
    
    pheno_month <- pheno_year %>% filter(Month_Num == month_num)
    
    species_productivity <- pheno_month %>%
      group_by(Species) %>%
      summarise(
        mean_fruit_per_tree = mean(`No. fruit on tree`),
        mean_seed_per_tree = mean(`No. seed pods on tree`),
        .groups = "drop"
      )
    
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
    
    veg_month$total_trees <- rowSums(veg_month[, species_cols])
    
    # FRUIT
    if (sum(veg_month$estimated_fruit > 0) >= 20) {
      
      veg_fruit <- veg_month %>% filter(estimated_fruit > 0)
      
      gam_fruit <- gam(
        log(estimated_fruit + 1) ~ s(UTM_X, UTM_Y, k = 12) +
          s(`Canopy.Cover(%)`, k = 5) +
          s(total_trees, k = 5),
        data = veg_fruit, method = "REML"
      )
      
      fruit_raster <- predict_to_raster(
        gam_fruit,
        raster_template,
        extra_vars = list(
          `Canopy.Cover(%)` = canopy_values,
          total_trees = mean(veg_month$total_trees)
        )
      )
      
      # Back-transform from log
      values(fruit_raster) <- exp(values(fruit_raster)) - 1
      
      dir.create(paste0("raster_outputs/", year), showWarnings = FALSE)
      writeRaster(fruit_raster, 
                  paste0("raster_outputs/", year, "/fruit_", 
                         tolower(month_name), ".tif"),
                  overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      
      rm(gam_fruit, fruit_raster)
      gc()
    }
    
    # SEED
    if (sum(veg_month$estimated_seeds > 0) >= 20) {
      
      veg_seed <- veg_month %>% filter(estimated_seeds > 0)
      
      gam_seed <- gam(
        log(estimated_seeds + 1) ~ s(UTM_X, UTM_Y, k = 12) +
          s(`Canopy.Cover(%)`, k = 5) +
          s(total_trees, k = 5),
        data = veg_seed, method = "REML"
      )
      
      seed_raster <- predict_to_raster(
        gam_seed,
        raster_template,
        extra_vars = list(
          `Canopy.Cover(%)` = canopy_values,
          total_trees = mean(veg_month$total_trees)
        )
      )
      
      # Back-transform from log
      values(seed_raster) <- exp(values(seed_raster)) - 1
      
      writeRaster(seed_raster, 
                  paste0("raster_outputs/", year, "/seed_", 
                         tolower(month_name), ".tif"),
                  overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      
      rm(gam_seed, seed_raster)
      gc()
    }
  }
}

cat("\n=== COMPLETE! ===\n")

# Summary
fruit_rasters <- list.files("raster_outputs", pattern = "fruit", recursive = TRUE)
seed_rasters <- list.files("raster_outputs", pattern = "seed", recursive = TRUE)

cat("Fruit rasters:", length(fruit_rasters), "\n")
cat("Seed rasters:", length(seed_rasters), "\n")
