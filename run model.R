

library(inlabru)
library(INLA)
library(dplyr)
library(amt)
library(raster)
library(sf)

#-- Load rasters

# Get DEM and Canopy height
dem <- raster("movement_analysis/dem.tif")

ext <- extent(dem)
bbox_utm <- st_bbox(c(xmin = ext@xmin, xmax = ext@xmax, 
                      ymin = ext@ymin, ymax = ext@ymax), 
                    crs = st_crs(32735))
bbox_latlon <- st_transform(st_as_sfc(bbox_utm), 4326)
st_bbox(bbox_latlon)

# Obtained .tif from GEE and manually imported into R
canopy_height <- raster("meta_tree_height.tif")
plot(canopy_height)

canopy_height_raster <- crop(canopy_height, dem)
plot(canopy_height_raster)

# Load canopy cover raster (needed for viewshed permeability)
cc_raster <- raster("canopy_cover.tif")

#-- Prepare GPS data

if (!("x_utm" %in% names(vervet_clean))) {
  vervet_sf <- st_as_sf(vervet_clean, coords = c("longitude", "latitude"), crs = 4326)
  vervet_utm <- st_transform(vervet_sf, crs = 32735)
  coords <- st_coordinates(vervet_utm)
  vervet_clean$x_utm <- coords[, 1]
  vervet_clean$y_utm <- coords[, 2]
}

if (!("x_utm" %in% names(baboon_clean))) {
  baboon_sf <- st_as_sf(baboon_clean, coords = c("longitude", "latitude"), crs = 4326)
  baboon_utm <- st_transform(baboon_sf, crs = 32735)
  coords <- st_coordinates(baboon_utm)
  baboon_clean$x_utm <- coords[, 1]
  baboon_clean$y_utm <- coords[, 2]
}

if ("New_Timestamp" %in% names(vervet_clean)) {
  vervet_clean$timestamp <- vervet_clean$New_Timestamp
  baboon_clean$timestamp <- baboon_clean$New_Timestamp
}

gps_dates <- unique(as.Date(vervet_clean$timestamp, tz = "UTC"))

cat("GPS data prepared\n")
cat("  Vervet:", nrow(vervet_clean), "locations\n")
cat("  Unique GPS dates:", length(gps_dates), "\n\n")

#-- Load and parse scan data

scans <- read.csv("scans_com.csv", stringsAsFactors = FALSE, quote = "\"")

scan_date_col <- names(scans)[grepl("Date", names(scans), ignore.case = TRUE)][1]
scan_time_col <- names(scans)[grepl("Time", names(scans), ignore.case = TRUE)][1]
scan_height_col <- names(scans)[grepl("Height", names(scans), ignore.case = TRUE)][1]

#-- Height parsing - handle "UN" as NA

clean_height <- function(h) {
  h <- trimws(as.character(h))
  
  # Handle "UN" as NA (unknown height)
  if (toupper(h) == "UN" || h == "" || is.na(h)) {
    return(NA)
  }
  
  # Handle >10
  if (grepl("^>\\s*10", h, ignore.case = TRUE)) {
    return(11)
  }
  
  # Remove non-numeric characters
  h_clean <- gsub("[^0-9.]", "", h)
  height <- as.numeric(h_clean)
  
  # Code >10 as 11
  if (!is.na(height) && height > 10) {
    return(11)
  }
  
  return(height)
}

scans$height_clean <- sapply(scans[[scan_height_col]], clean_height)

cat("  Total rows:", nrow(scans), "\n")
cat("  Valid heights:", sum(!is.na(scans$height_clean)), "\n")
cat("  NA/UN heights:", sum(is.na(scans$height_clean)), "\n\n")

#-- Date parsing (GPS-informed)

parse_date_gps_informed <- function(date_str, gps_dates_ref) {
  date_str <- trimws(as.character(date_str))
  
  # Numeric dates
  if (grepl("^\\d+$", date_str)) {
    r_numeric <- as.numeric(date_str)
    date <- as.Date(r_numeric, origin = "1970-01-01")
    year <- as.integer(format(date, "%Y"))
    if (!is.na(date) && year >= 2014 && year <= 2019) {
      return(date)
    }
    return(NA)
  }
  
  # Text dates
  if (grepl("/", date_str)) {
    parts <- strsplit(date_str, "/")[[1]]
    if (length(parts) != 3) return(NA)
    
    part1 <- as.integer(parts[1])
    part2 <- as.integer(parts[2])
    year <- as.integer(parts[3])
    
    # Unambiguous
    if (part1 > 12) {
      return(as.Date(paste(year, part2, part1, sep = "-")))
    }
    if (part2 > 12) {
      return(as.Date(paste(year, part1, part2, sep = "-")))
    }
    
    # Ambiguous - use GPS reference
    uk_date <- tryCatch(as.Date(paste(year, part2, part1, sep = "-")), error = function(e) NA)
    us_date <- tryCatch(as.Date(paste(year, part1, part2, sep = "-")), error = function(e) NA)
    
    uk_in_gps <- !is.na(uk_date) && (uk_date %in% gps_dates_ref)
    us_in_gps <- !is.na(us_date) && (us_date %in% gps_dates_ref)
    
    if (uk_in_gps && !us_in_gps) return(uk_date)
    if (us_in_gps && !uk_in_gps) return(us_date)
    if (!is.na(uk_date)) return(uk_date)  # Default UK
    return(us_date)
  }
  
  return(NA)
}

scans$date_clean <- sapply(scans[[scan_date_col]], parse_date_gps_informed, gps_dates_ref = gps_dates)
scans$date_clean <- as.Date(scans$date_clean, origin = "1970-01-01")


#-- Time
parse_time_robust <- function(time_str) {
  time_str <- trimws(as.character(time_str))
  if (grepl("^\\d{1,2}:\\d{2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:%02d", as.integer(parts[1]), as.integer(parts[2]), as.integer(parts[3])))
  }
  if (grepl("^\\d{1,2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:00", as.integer(parts[1]), as.integer(parts[2])))
  }
  return(NA)
}

scans$time_clean <- sapply(scans[[scan_time_col]], parse_time_robust)

#-- Create datetime
scans$datetime <- as.POSIXct(paste(scans$date_clean, scans$time_clean), 
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


# First filter: keep rows with valid datetime
scans_valid <- scans %>%
  filter(!is.na(datetime))

# Group by scan event (date + time)
scan_samples <- scans_valid %>%
  group_by(date_clean, datetime) %>%
  summarize(
    n_individuals = n(),
    n_valid_heights = sum(!is.na(height_clean)),
    max_height_m = ifelse(n_valid_heights > 0, 
                          max(height_clean, na.rm = TRUE), 
                          NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(max_height_m)) %>%  # Drop scan samples where ALL heights are NA
  arrange(datetime)



# SSF steps
ver_track <- vervet_clean %>%
  make_track(longitude, latitude, timestamp, crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)

ver_steps <- ver_track %>% steps()

# Filter out overnight steps and long gaps (>2 hours)
cat("  Steps before filtering:", nrow(ver_steps), "\n")

ver_steps <- ver_steps %>%
  mutate(
    # Calculate time difference in hours
    dt_hours = as.numeric(difftime(t2_, t1_, units = "hours")),
    
    # Check if same date
    date1 = as.Date(t1_, tz = "UTC"),
    date2 = as.Date(t2_, tz = "UTC"),
    same_date = (date1 == date2)
  ) %>%
  # Keep only steps within same day AND <= 2 hours
  filter(same_date == TRUE, dt_hours <= 2.0) %>%
  dplyr::select(-date1, -date2, -same_date)  # Clean up temp columns

cat("  Steps after filtering overnight/long gaps:", nrow(ver_steps), "\n")

ver_ssf_data <- ver_steps %>% random_steps(n_control = 10)

cat("  Total SSF rows:", nrow(ver_ssf_data), "\n")
cat("  Observed:", sum(ver_ssf_data$case_), "\n")
cat("  Random:", sum(!ver_ssf_data$case_), "\n\n")


#-- Match sentinel heights

# Max height per scan window
strata_info <- ver_ssf_data %>%
  filter(case_ == TRUE) %>%
  dplyr::select(step_id_, t2_)

sentinel_heights <- numeric(nrow(strata_info))
height_source <- character(nrow(strata_info))
time_diff_hours <- numeric(nrow(strata_info))

for (i in 1:nrow(strata_info)) {
  step_time <- strata_info$t2_[i]
  step_date <- as.Date(step_time, tz = "UTC")
  
  same_day_scans <- scan_samples %>% filter(date_clean == step_date)
  
  if (nrow(same_day_scans) > 0) {
    time_diffs <- abs(difftime(same_day_scans$datetime, step_time, units = "hours"))
    nearest_idx <- which.min(time_diffs)
    
    sentinel_heights[i] <- same_day_scans$max_height_m[nearest_idx]
    time_diff_hours[i] <- as.numeric(time_diffs[nearest_idx])
    height_source[i] <- "scan_same_day"
  } else {
    sentinel_heights[i] <- NA
    time_diff_hours[i] <- NA
    height_source[i] <- "no_scan_on_date"
  }
}

strata_info$sentinel_height_m <- sentinel_heights
strata_info$height_source <- height_source
strata_info$time_diff_hours <- time_diff_hours

# Merge back to all steps
ver_ssf_data <- ver_ssf_data %>%
  left_join(strata_info %>% dplyr::select(step_id_, sentinel_height_m, height_source, time_diff_hours), 
            by = "step_id_")

print(table(height_source))

# Filter heights with na
# Identify strata with valid heights
valid_strata <- strata_info %>%
  filter(!is.na(sentinel_height_m)) %>%
  pull(step_id_)

# Keep only steps from valid strata
ver_ssf_data <- ver_ssf_data %>%
  filter(step_id_ %in% valid_strata)

cat("  Steps after filtering:", nrow(ver_ssf_data), "\n")
cat("  Strata retained:", length(unique(ver_ssf_data$step_id_)), "\n")

#-- Viewshed clearance function
calculate_viewshed_clearance <- function(observer_x, observer_y, observer_height_agl,
                                         target_x, target_y, 
                                         dem, canopy_height_raster,
                                         target_height_agl = 1.5) {
  
  # Extract ground elevations
  obs_ground <- raster::extract(dem, matrix(c(observer_x, observer_y), ncol = 2))
  tgt_ground <- raster::extract(dem, matrix(c(target_x, target_y), ncol = 2))
  
  if (is.na(obs_ground) || is.na(tgt_ground)) {
    return(list(min_clearance = NA, distance_m = NA, n_obstructions = NA))
  }
  
  # Absolute elevations of observer and target
  observer_elev <- obs_ground + observer_height_agl
  target_elev <- tgt_ground + target_height_agl
  distance <- sqrt((target_x - observer_x)^2 + (target_y - observer_y)^2)
  
  # Sample points along line of sight
  n_points <- max(10, ceiling(distance / 10))
  x_points <- seq(observer_x, target_x, length.out = n_points)
  y_points <- seq(observer_y, target_y, length.out = n_points)
  
  # Extract terrain and vegetation data
  terrain_elevs <- raster::extract(dem, cbind(x_points, y_points))
  canopy_heights <- raster::extract(canopy_height_raster, cbind(x_points, y_points))
  canopy_heights[is.na(canopy_heights)] <- 0
  
  # CRITICAL: Use FULL canopy heights - never reduce them!
  # Obstacle top elevation = terrain + FULL vegetation height
  obstacle_elevs <- terrain_elevs + canopy_heights
  
  # Calculate line-of-sight elevation at each point
  fractions <- seq(0, 1, length.out = n_points)
  los_elevs <- observer_elev + fractions * (target_elev - observer_elev)
  
  # Geometric clearance: vertical distance from LOS to obstacle tops
  # Positive = LOS passes ABOVE obstacles (clear)
  # Negative = LOS passes BELOW obstacles (obstructed)
  clearances <- los_elevs - obstacle_elevs
  
  # Exclude endpoints (don't count observer/target as obstacles)
  if (length(clearances) > 2) {
    clearances <- clearances[2:(length(clearances)-1)]
  }
  
  # Minimum clearance is the most restrictive point (bottleneck)
  min_clearance <- min(clearances, na.rm = TRUE)
  n_obstructions <- sum(clearances < 0, na.rm = TRUE)
  
  return(list(
    min_clearance = min_clearance, 
    distance_m = distance, 
    n_obstructions = n_obstructions
  ))
}


# Conceptual framework:
# - viewshed_clearance: 3D geometric sight-line obstruction
# - canopy_cover_pct: permeability/detectability through vegetation
# - horizontal_visibility: ambient ground-level visibility

#-- Calculate covariates
ver_ssf_data$dist_to_baboon <- NA
ver_ssf_data$viewshed_clearance_m <- NA
ver_ssf_data$n_obstructions <- NA

start_time <- Sys.time()
n_processed <- 0

for (i in 1:nrow(ver_ssf_data)) {
  step_x <- ver_ssf_data$x2_[i]
  step_y <- ver_ssf_data$y2_[i]
  step_time <- ver_ssf_data$t2_[i]
  sentinel_height <- ver_ssf_data$sentinel_height_m[i]
  
  time_diffs <- abs(difftime(baboon_clean$timestamp, step_time, units = "mins"))
  nearest_bab_idx <- which.min(time_diffs)
  
  if (length(nearest_bab_idx) == 0 || time_diffs[nearest_bab_idx] > 10) next
  
  baboon_x <- baboon_clean$x_utm[nearest_bab_idx]
  baboon_y <- baboon_clean$y_utm[nearest_bab_idx]
  
  dist <- sqrt((step_x - baboon_x)^2 + (step_y - baboon_y)^2)
  ver_ssf_data$dist_to_baboon[i] <- dist
  
  # Calculate viewshed clearance WITHOUT canopy cover adjustment
  vs <- calculate_viewshed_clearance(step_x, step_y, sentinel_height,
                                     baboon_x, baboon_y, 
                                     dem, canopy_height_raster)
  
  ver_ssf_data$viewshed_clearance_m[i] <- vs$min_clearance
  ver_ssf_data$n_obstructions[i] <- vs$n_obstructions
  
  n_processed <- n_processed + 1
  
  if (n_processed %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- n_processed / elapsed
    remaining <- (nrow(ver_ssf_data) - n_processed) / rate
    cat(sprintf("  %d/%d (%.1f%%) - %.1f steps/sec - ETA: %.1f min\r",
                n_processed, nrow(ver_ssf_data), 100*n_processed/nrow(ver_ssf_data), 
                rate, remaining/60))
  }
}


#--
#-- Extract additional covariates
#--

# Extract from horizontal visibility raster
hv_raster <- raster("horizontal_visibility.tif")
coords_matrix <- cbind(ver_ssf_data$x2_, ver_ssf_data$y2_)
ver_ssf_data$horizontal_visibility_m <- raster::extract(hv_raster, coords_matrix)
ver_ssf_data$canopy_cover_pct <- raster::extract(cc_raster, coords_matrix)

write.csv(ver_ssf_data, "ver_ssf_data_complete.csv", row.names = FALSE)
str(ver_ssf_data)

#--
#-- Prepare for modeling
#--

ssf_complete <- ver_ssf_data %>%
  filter(!is.na(dist_to_baboon), !is.na(viewshed_clearance_m))

# Scale covariates
ssf_complete$dist_scaled <- scale(ssf_complete$dist_to_baboon)[,1]
ssf_complete$clearance_scaled <- scale(ssf_complete$viewshed_clearance_m)[,1]
ssf_complete$vis_scaled <- scale(ssf_complete$horizontal_visibility_m)[,1]

# Calculate weights
ssf_complete <- ssf_complete %>%
  group_by(step_id_) %>%
  mutate(
    n_steps_in_stratum = n(),
    weight = 1e6 / n_steps_in_stratum
  ) %>%
  ungroup()

str(ssf_complete)

# Convert to integer for INLA
ssf_complete$case_numeric <- as.integer(ssf_complete$case_)
ssf_complete$stratum_id <- as.factor(ssf_complete$step_id_)

# Scale canopy cover for interaction
ssf_complete$canopy_scaled <- scale(ssf_complete$canopy_cover_pct)[,1]

# Create all distance interactions
ssf_complete$dist_x_clear <- ssf_complete$dist_scaled * ssf_complete$clearance_scaled
ssf_complete$dist_x_canopy <- ssf_complete$dist_scaled * ssf_complete$canopy_scaled
ssf_complete$dist_x_vis <- ssf_complete$dist_scaled * ssf_complete$vis_scaled

#-- Fit model
global_model <- inla(
  case_numeric ~ -1 + 
    # Main effects (all fixed in model)
    dist_scaled + 
    clearance_scaled +
    vis_scaled + 
    
    # Two-way interactions (test detection-mediated avoidance)
    dist_x_clear +
    
    # Stratum random effects
    f(stratum_id, model = "iid", hyper = list(prec = list(initial = -6, fixed = TRUE))),
  
  family = "poisson",
  data = ssf_complete,
  E = weight,
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  verbose = FALSE
)

summary(global_model)

#--
#-- Display results
#--

fixed <- global_model$summary.fixed

for (i in 1:nrow(fixed)) {
  param <- rownames(fixed)[i]
  beta <- fixed[i, "mean"]
  ci_low <- fixed[i, "0.025quant"]
  ci_high <- fixed[i, "0.975quant"]
  
  # Denote if Credible interval includes 0
  sig <- if (ci_low > 0 | ci_high < 0) " *" else ""
  
  cat(sprintf("%-25s β = %7.3f  (95%% CI: [%7.3f, %7.3f])%s\n",
              param, beta, ci_low, ci_high, sig))
}

cat(sprintf("DIC:  %.1f\n", global_model$dic$dic))
cat(sprintf("WAIC: %.1f\n", global_model$waic$waic))
cat(sprintf("Effective parameters: %.1f\n\n", global_model$dic$p.eff))

#--
# Save results
#--

results_summary <- list(
  fixed_effects = as.data.frame(fixed),
  dic = global_model$dic$dic,
  waic = global_model$waic$waic,
  n_parameters = global_model$dic$p.eff
)

saveRDS(results_summary, "model_results_summary.rds")
