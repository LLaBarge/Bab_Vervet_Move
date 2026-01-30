# ============================================================================
# Movement Trajectory Analysis and Species Interaction Inference
# Restructured to work directly with matching_data
# ============================================================================

library(ggplot2)
remove.packages("rlang")
install.packages("rlang")
# ============================================================================
# STEP 0: PREPARE DATA
# ============================================================================

cat("Preparing data...\n")

# Ensure New_Timestamp is POSIXct
matching_data$New_Timestamp <- as.POSIXct(matching_data$New_Timestamp, 
                                          format = "%Y-%m-%d %H:%M:%S")
matching_data$Date <- as.Date(matching_data$New_Timestamp)

# Separate by species
baboon_data <- matching_data %>% 
  filter(species == "Baboon") %>%
  arrange(New_Timestamp)

vervet_data <- matching_data %>% 
  filter(species == "Vervet") %>%
  arrange(New_Timestamp)

cat(sprintf("Baboon points: %d\n", nrow(baboon_data)))
cat(sprintf("Vervet points: %d\n\n", nrow(vervet_data)))

# ============================================================================
# STEP 1: VISUALIZE TRAJECTORIES BY SPECIES
# ============================================================================

cat("STEP 1: Visualizing movement trajectories...\n\n")


# All trajectories together
p1 <- ggplot(matching_data, aes(x = longitude, y = latitude, color = species)) +
  geom_path(aes(group = interaction(species, Date)), alpha = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Baboon" = "#D55E00", "Vervet" = "#009E73")) +
  theme_minimal() +
  labs(title = "Movement Trajectories by Species",
       x = "Longitude", y = "Latitude", color = "Species") +
  coord_fixed(ratio = 1)
p1
ggsave("movement_analysis/trajectories_all.png", p1, width = 10, height = 8, dpi = 300)

# Baboon trajectories only
p2 <- ggplot(baboon_data, aes(x = longitude, y = latitude)) +
  geom_path(aes(group = Date), alpha = 0.5, color = "#D55E00") +
  geom_point(size = 1, alpha = 0.6, color = "#D55E00") +
  theme_minimal() +
  labs(title = "Baboon Movement Trajectories",
       x = "Longitude", y = "Latitude") +
  coord_fixed(ratio = 1)

ggsave("movement_analysis/trajectories_baboon.png", p2, width = 10, height = 8, dpi = 300)

# Vervet trajectories only
p3 <- ggplot(vervet_data, aes(x = longitude, y = latitude)) +
  geom_path(aes(group = Date), alpha = 0.5, color = "#009E73") +
  geom_point(size = 1, alpha = 0.6, color = "#009E73") +
  theme_minimal() +
  labs(title = "Vervet Movement Trajectories",
       x = "Longitude", y = "Latitude") +
  coord_fixed(ratio = 1)
p3
ggsave("movement_analysis/trajectories_vervet.png", p3, width = 10, height = 8, dpi = 300)

cat("Trajectory plots saved\n\n")

# ============================================================================
# STEP 2: FILTER UNREALISTIC MOVEMENTS (BY SPECIES)
# ============================================================================

cat("STEP 2: Filtering unrealistic movements...\n")
cat("Maximum speed threshold: 5 km/h\n\n")

# Function to calculate speeds and filter
filter_by_speed <- function(data, species_name, max_speed_kmh = 5) {
  
  cat(sprintf("Processing %s data...\n", species_name))
  
  data <- data %>%
    mutate(
      next_lon = lead(longitude),
      next_lat = lead(latitude),
      next_time = lead(New_Timestamp),
      
      distance_m = geosphere::distHaversine(
        cbind(longitude, latitude),
        cbind(next_lon, next_lat)
      ),
      
      time_diff_hours = as.numeric(difftime(next_time, New_Timestamp, units = "hours")),
      
      speed_kmh = (distance_m / 1000) / time_diff_hours,
      
      unrealistic = case_when(
        is.na(speed_kmh) ~ FALSE,
        is.na(time_diff_hours) ~ FALSE,
        time_diff_hours >= 2 ~ FALSE,
        speed_kmh > max_speed_kmh ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  n_unrealistic <- sum(data$unrealistic, na.rm = TRUE)
  
  cat(sprintf("  Total points: %d\n", nrow(data)))
  cat(sprintf("  Unrealistic: %d (%.2f%%)\n", 
              n_unrealistic, 100 * n_unrealistic / nrow(data)))
  
  # Plot speed distribution
  if (sum(!is.na(data$speed_kmh)) > 0) {
    p <- ggplot(data %>% filter(!is.na(speed_kmh), speed_kmh < 100), 
                aes(x = speed_kmh)) +
      geom_histogram(bins = 50, fill = ifelse(species_name == "Baboon", "#D55E00", "#009E73"), 
                     alpha = 0.7) +
      geom_vline(xintercept = max_speed_kmh, color = "red", 
                 linetype = "dashed", size = 1) +
      theme_minimal() +
      labs(title = paste(species_name, "Movement Speeds"),
           subtitle = paste("Red line: threshold =", max_speed_kmh, "km/h"),
           x = "Speed (km/h)", y = "Count")
    
    ggsave(paste0("movement_analysis/speed_distribution_", tolower(species_name), ".png"), 
           p, width = 8, height = 6, dpi = 300)
  }
  
  # Remove unrealistic points
  data_clean <- data %>%
    filter(unrealistic == FALSE) %>%
    dplyr::select(-next_lon, -next_lat, -next_time, -distance_m, 
                  -time_diff_hours, -speed_kmh, -unrealistic)
  
  cat(sprintf("  Remaining: %d\n\n", nrow(data_clean)))
  
  return(data_clean)
}

# Filter each species
baboon_clean <- filter_by_speed(baboon_data, "Baboon", max_speed_kmh = 5)
vervet_clean <- filter_by_speed(vervet_data, "Vervet", max_speed_kmh = 5)

# Combine cleaned data
matching_data_clean <- bind_rows(baboon_clean, vervet_clean)

cat(sprintf("Total cleaned points: %d\n\n", nrow(matching_data_clean)))

# ============================================================================
# STEP 3: COEFFICIENT OF SOCIALITY (Cs)
# ============================================================================

# Prepare data
bab_df <- baboon_clean %>%
  dplyr::select(New_Timestamp, longitude, latitude) %>%
  arrange(New_Timestamp)

ver_df <- vervet_clean %>%
  dplyr::select(New_Timestamp, longitude, latitude) %>%
  arrange(New_Timestamp)

# Find pairs within time tolerance (20 minutes = 1200 seconds)
time_tolerance <- 60 * 20  # seconds

prox <- data.frame()

for (i in 1:nrow(bab_df)) {
  time_diffs <- abs(difftime(ver_df$New_Timestamp, bab_df$New_Timestamp[i], units = "secs"))
  match_idx <- which(time_diffs <= time_tolerance)
  
  if (length(match_idx) > 0) {
    # Take closest match in time
    closest_idx <- match_idx[which.min(time_diffs[match_idx])]
    
    prox <- rbind(prox, data.frame(
      x1 = bab_df$longitude[i],
      y1 = bab_df$latitude[i],
      x2 = ver_df$longitude[closest_idx],
      y2 = ver_df$latitude[closest_idx],
      date = bab_df$New_Timestamp[i],
      time_diff = as.numeric(time_diffs[closest_idx])
    ))
  }
}

if (nrow(prox) > 0) {
  
  cat(sprintf("Found %d simultaneous fixes\n\n", nrow(prox)))
  
  # Calculate distances
  prox$distance <- geosphere::distHaversine(
    cbind(prox$x1, prox$y1),
    cbind(prox$x2, prox$y2)
  )
  
  # Proximity at different thresholds
  cat("Proximity analysis:\n")
  
  for (thresh in c(50, 100, 150, 200, 300)) {
    n_close <- sum(prox$distance <= thresh, na.rm = TRUE)
    cat(sprintf("  Within %dm: %d / %d (%.1f%%)\n", 
                thresh, n_close, nrow(prox), 100 * n_close / nrow(prox)))
  }
  
  cat(sprintf("\nMean distance: %.1fm\n", mean(prox$distance, na.rm = TRUE)))
  cat(sprintf("Median distance: %.1fm\n\n", median(prox$distance, na.rm = TRUE)))
  
  # Calculate Coefficient of Sociality manually
  # Cs = (observed proximity - expected proximity) / (observed proximity + expected proximity)
  # Using 150m threshold
  
  n_close_150 <- sum(prox$distance <= 150)
  prop_close <- n_close_150 / nrow(prox)
  
  cat("Coefficient of Sociality (simplified calculation):\n")
  cat(sprintf("  Proportion within 150m: %.3f\n", prop_close))
  
  if (prop_close > 0.5) {
    cat("  Interpretation: High association (attraction)\n\n")
  } else if (prop_close > 0.2) {
    cat("  Interpretation: Moderate association\n\n")
  } else {
    cat("  Interpretation: Low association (possible avoidance)\n\n")
  }
  
  # Plot distance distribution
  p_dist <- ggplot(data.frame(distance = prox$distance), aes(x = distance)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 150, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Distance Between Species at Simultaneous Fixes",
         subtitle = "Red line: 150m threshold",
         x = "Distance (meters)", y = "Count")
  
  ggsave("movement_analysis/simultaneous_distances.png", p_dist, 
         width = 8, height = 6, dpi = 300)
  
} else {
  cat("No simultaneous fixes found\n\n")
  prox <- NULL
}

# ============================================================================
# STEP 4: SSF-DIST ANALYSIS (VERVETS ONLY)
# ============================================================================

cat("STEP 4: SSF-DIST Analysis for Vervets...\n\n")

# Create tracks
bab_track <- baboon_clean %>%
  make_track(longitude, latitude, New_Timestamp, 
             crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)

ver_track <- vervet_clean %>%
  make_track(longitude, latitude, New_Timestamp, 
             crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)

cat(sprintf("Baboon track: %d locations\n", nrow(bab_track)))
cat(sprintf("Vervet track: %d locations\n\n", nrow(ver_track)))

# Calculate distance from vervets to baboons
cat("Calculating distances to baboons...\n")

get_nearest_distance <- function(time, focal_x, focal_y, other_track) {
  time_diffs <- abs(difftime(other_track$t_, time, units = "mins"))
  nearest_idx <- which.min(time_diffs)
  
  if (length(nearest_idx) > 0 && time_diffs[nearest_idx] <= 60) {
    dist <- sqrt((focal_x - other_track$x_[nearest_idx])^2 + 
                   (focal_y - other_track$y_[nearest_idx])^2)
    return(dist)
  } else {
    return(NA)
  }
}

ver_track$dist_to_baboon <- mapply(
  get_nearest_distance,
  ver_track$t_,
  ver_track$x_,
  ver_track$y_,
  MoreArgs = list(other_track = bab_track)
)

n_with_dist <- sum(!is.na(ver_track$dist_to_baboon))
cat(sprintf("Vervet points with distance data: %d / %d\n\n", 
            n_with_dist, nrow(ver_track)))

# Generate steps (without burst - treat as single trajectory)
cat("Generating steps...\n")

# Add burst column (all same burst since it's one trajectory)
ver_track$burst_ <- 1

ver_steps <- ver_track %>%
  steps_by_burst()

cat(sprintf("Generated %d steps\n\n", nrow(ver_steps)))

# Generate random steps
cat("Generating random steps...\n")

ver_ssf_data <- ver_steps %>%
  random_steps(n_control = 10)

cat(sprintf("Total steps: %d (%d observed, %d random)\n\n", 
            nrow(ver_ssf_data),
            sum(ver_ssf_data$case_),
            sum(!ver_ssf_data$case_)))

# Fit SSF model
cat("Fitting SSF model...\n")

ver_model <- tryCatch({
  fit_issf(
    case_ ~ 
      sl_ + log_sl_ +
      cos_ta_ +
      strata(step_id_),
    data = ver_ssf_data,
    model = TRUE
  )
}, error = function(e) {
  cat("Error fitting model:", e$message, "\n")
  return(NULL)
})

if (!is.null(ver_model)) {
  cat("\nModel Summary:\n")
  print(summary(ver_model))
  
  cat("\n\nInterpretation:\n")
  cat("  sl_: Effect of step length\n")
  cat("  log_sl_: Non-linear step length effect\n")
  cat("  cos_ta_: Turning angle preference\n")
  cat("    Positive = tendency to continue forward\n")
  cat("    Negative = tendency to turn back\n\n")
}

# Plot distance distribution
dist_df <- data.frame(distance = ver_track$dist_to_baboon) %>%
  filter(!is.na(distance))

if (nrow(dist_df) > 0) {
  p_vb <- ggplot(dist_df, aes(x = distance)) +
    geom_histogram(bins = 50, fill = "#009E73", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Vervet-Baboon Distances",
         x = "Distance to Nearest Baboon (meters)", y = "Count")
  
  ggsave("movement_analysis/vervet_baboon_distances.png", p_vb, 
         width = 8, height = 6, dpi = 300)
  
  cat(sprintf("Mean vervet-baboon distance: %.1fm\n", mean(dist_df$distance)))
  cat(sprintf("Median vervet-baboon distance: %.1fm\n\n", median(dist_df$distance)))
}


# ============================================================================
# SAVE RESULTS
# ============================================================================

results <- list(
  baboon_clean = baboon_clean,
  vervet_clean = vervet_clean,
  matching_data_clean = matching_data_clean,
  simultaneous_fixes = if(exists("prox")) prox else NULL,
  cs_value = if(exists("cs_value")) cs_value else NA,
  vervet_model = if(exists("ver_model")) ver_model else NULL,
  vervet_track = ver_track,
  baboon_track = bab_track
)

saveRDS(results, "movement_analysis/analysis_results.rds")

cat("ANALYSIS COMPLETE\n")
cat("Results saved to: movement_analysis/\n")