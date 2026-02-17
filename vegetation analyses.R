#===============================================================================
# vegetation data
#===============================================================================
setwd("C:/Users/lrlab/PPP scan points")

library(readxl)
library(dplyr)
library(sf)
library(terra)
library(mgcv)

# Load data
veg_plots  <- read_excel("vegplot_phenology.xlsx", sheet = 1)
phenology  <- read_excel("vegplot_phenology.xlsx", sheet = "phenology")

# Replace NA with 0
veg_plots[is.na(veg_plots)] <- 0
phenology[is.na(phenology)]  <- 0

# Phenology
phenology$Month <- tools::toTitleCase(tolower(trimws(phenology$Month)))
month_mapping <- c(
  "January"=1,"February"=2,"March"=3,"April"=4,
  "May"=5,"June"=6,"July"=7,"August"=8,
  "September"=9,"October"=10,"November"=11,"December"=12
)
phenology$Month_Num <- month_mapping[phenology$Month]
phenology_filtered  <- phenology %>% filter(Year >= 2014 & Year <= 2018)
phenology_filtered$`No. fruit on tree`    <- as.numeric(phenology_filtered$`No. fruit on tree`)
phenology_filtered$`No. seed pods on tree`<- as.numeric(phenology_filtered$`No. seed pods on tree`)
phenology_filtered[is.na(phenology_filtered)] <- 0

# Vegetation plot species columns and summaries
species_cols           <- names(veg_plots)[grepl("^[A-Z]\\.", names(veg_plots))]
veg_plots$total_trees  <- rowSums(veg_plots[, species_cols])
veg_plots$species_richness <- rowSums(veg_plots[, species_cols] > 0)

# UTM coordinates
veg_sf  <- st_as_sf(veg_plots, coords = c("Long", "Lat"), crs = 4326)
veg_utm <- st_transform(veg_sf, crs = 32736)
coords  <- st_coordinates(veg_utm)
veg_plots$UTM_X <- coords[, 1]
veg_plots$UTM_Y <- coords[, 2]

veg_plots$CH <- veg_plots$`canopy.height(m)`
veg_plots$CC <- veg_plots$`Canopy.Cover(%)`

# rast template
RESOLUTION <- 100   # metres

x_range  <- range(veg_plots$UTM_X)
y_range  <- range(veg_plots$UTM_Y)
x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05

raster_template <- rast(
  xmin       = x_range[1] - x_buffer,
  xmax       = x_range[2] + x_buffer,
  ymin       = y_range[1] - y_buffer,
  ymax       = y_range[2] + y_buffer,
  resolution = RESOLUTION,
  crs        = "EPSG:32736"
)

cat("Resolution:", RESOLUTION, "m\n")
cat("Raster size:", nrow(raster_template), "x", ncol(raster_template), "\n")
cat("Total cells:", ncell(raster_template), "\n\n")

# predict to rasters
predict_to_raster <- function(model, template, extra_vars = NULL) {
  
  coords_matrix <- crds(template, df = FALSE)
  chunk_size    <- 10000   
  n_cells       <- nrow(coords_matrix)
  n_chunks      <- ceiling(n_cells / chunk_size)
  predictions   <- numeric(n_cells)
  
  for (i in 1:n_chunks) {
    start_idx  <- (i - 1) * chunk_size + 1
    end_idx    <- min(i * chunk_size, n_cells)
    
    chunk_data <- data.frame(
      UTM_X = coords_matrix[start_idx:end_idx, 1],
      UTM_Y = coords_matrix[start_idx:end_idx, 2]
    )
    
    if (!is.null(extra_vars)) {
      for (var_name in names(extra_vars)) {
        chunk_data[[var_name]] <- extra_vars[[var_name]]
      }
    }
    
    predictions[start_idx:end_idx] <- predict(model, newdata = chunk_data,
                                              type = "response")
    
    if (i %% 20 == 0) cat("  Chunk", i, "of", n_chunks, "\n")
  }
  
  values(template) <- predictions
  return(template)
}

# ============================================================================
# STATIC RASTERS: CANOPY COVER AND CANOPY HEIGHT
# ============================================================================

dir.create("raster_outputs", showWarnings = FALSE)

cat("Fitting canopy cover GAM...\n")
gam_canopy_cover <- gam(CC ~ s(UTM_X, UTM_Y, k = 25),
                        data = veg_plots, method = "REML")
gc()
cat("Predicting canopy cover raster...\n")
canopy_cover_raster <- predict_to_raster(gam_canopy_cover, raster_template)
writeRaster(canopy_cover_raster, "raster_outputs/canopy_cover.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("Fitting canopy height GAM...\n")
gam_canopy_height <- gam(CH ~ s(UTM_X, UTM_Y, k = 25),
                         data = veg_plots, method = "REML")

cat("Predicting canopy height raster...\n")
canopy_height_raster <- predict_to_raster(gam_canopy_height, raster_template)
writeRaster(canopy_height_raster, "raster_outputs/canopy_height.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

rm(gam_canopy_cover, gam_canopy_height)
gc()
cat("Static rasters done.\n\n")

# ============================================================================
# MONTHLY FRUIT AND SEED RASTERS
# ============================================================================

years <- sort(unique(phenology_filtered$Year))

for (year in years) {
  
  dir.create(paste0("raster_outputs/", year), showWarnings = FALSE)
  
  pheno_year     <- phenology_filtered %>% filter(Year == year)
  months_in_year <- sort(unique(pheno_year$Month_Num))
  
  for (month_num in months_in_year) {
    
    month_name <- names(month_mapping)[month_mapping == month_num]
    cat(month_name, "...\n")
    
    pheno_month <- pheno_year %>% filter(Month_Num == month_num)
    
    species_productivity <- pheno_month %>%
      group_by(Species) %>%
      summarise(
        mean_fruit_per_tree = mean(`No. fruit on tree`),
        mean_seed_per_tree  = mean(`No. seed pods on tree`),
        .groups = "drop"
      )
    
    veg_month <- veg_plots
    veg_month$estimated_fruit <- 0
    veg_month$estimated_seeds <- 0
    
    for (species in species_cols) {
      if (species %in% species_productivity$Species) {
        sp_data <- species_productivity %>% filter(Species == species)
        veg_month$estimated_fruit <- veg_month$estimated_fruit +
          veg_month[[species]] * sp_data$mean_fruit_per_tree
        veg_month$estimated_seeds <- veg_month$estimated_seeds +
          veg_month[[species]] * sp_data$mean_seed_per_tree
      }
    }
    
    veg_month$total_trees <- rowSums(veg_month[, species_cols])
    
    # ── FRUIT ────────────────────────────────────────────────
    if (sum(veg_month$estimated_fruit > 0) >= 20) {
      
      veg_fruit <- veg_month %>% filter(estimated_fruit > 0)
      
      gam_fruit <- gam(
        log(estimated_fruit + 1) ~ s(UTM_X, UTM_Y, k = 12) +
          s(CC, k = 5) +
          s(total_trees, k = 5),
        data = veg_fruit, method = "REML"
      )
      
      # FIX: total_trees passed inside extra_vars, not as direct argument
      fruit_raster <- predict_to_raster(
        gam_fruit,
        raster_template,
        extra_vars = list(
          CC          = mean(veg_month$CC),
          total_trees = mean(veg_month$total_trees)
        )
      )
      
      values(fruit_raster) <- exp(values(fruit_raster)) - 1
      
      writeRaster(fruit_raster,
                  paste0("raster_outputs/", year, "/fruit_",
                         tolower(month_name), ".tif"),
                  overwrite = TRUE, gdal = c("COMPRESS=LZW"))
      
      rm(gam_fruit, fruit_raster)
      gc()
    }
    
    # ── SEED ─────────────────────────────────────────────────
    if (sum(veg_month$estimated_seeds > 0) >= 20) {
      
      veg_seed <- veg_month %>% filter(estimated_seeds > 0)
      
      gam_seed <- gam(
        log(estimated_seeds + 1) ~ s(UTM_X, UTM_Y, k = 12) +
          s(CC, k = 5) +
          s(total_trees, k = 5),
        data = veg_seed, method = "REML"
      )
      
      # FIX: both CC and total_trees inside extra_vars
      seed_raster <- predict_to_raster(
        gam_seed,
        raster_template,
        extra_vars = list(
          CC          = mean(veg_month$CC),
          total_trees = mean(veg_month$total_trees)
        )
      )
      
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

# Summary
fruit_rasters <- list.files("raster_outputs", pattern = "fruit", recursive = TRUE)
seed_rasters  <- list.files("raster_outputs", pattern = "seed",  recursive = TRUE)
cat("\nFruit rasters:", length(fruit_rasters), "\n")
cat("Seed rasters: ", length(seed_rasters),  "\n")
