library(geosphere)
library(adehabitatLT)
devtools::install_github('thomasp85/gganimate', force = TRUE)


# Install and load required packages
required_packages <- c(
  "dplyr", "ggplot2", "sf", "lubridate","geosphere" "gifski",
  "move", "adehabitatLT", "wildlifeDI", "amt", "viridis", "adehabitatLT"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# PART 1: VISUALIZE MOVEMENT TRAJECTORIES BY DATE
# ============================================================================

visualize_trajectories <- function(data, output_dir = "outputs") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ensure New_Timestamp is POSIXct
  data$New_Timestamp <- as.POSIXct(data$New_Timestamp, format = "%Y-%m-%d %H:%M:%S")
  data$Date <- as.Date(data$New_Timestamp)
  
  # Filter out NA timestamps
  data <- data %>% filter(!is.na(New_Timestamp))
  
  # Static plot: All trajectories colored by species
  p1 <- ggplot(data, aes(x = longitude, y = latitude, color = species)) +
    geom_path(aes(group = interaction(species, Date)), alpha = 0.3) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = c("Baboon" = "#D55E00", "Vervet" = "#009E73")) +
    theme_minimal() +
    labs(title = "Movement Trajectories by Species",
         x = "Longitude", y = "Latitude",
         color = "Species") +
    coord_fixed(ratio = 1)
  
  ggsave(file.path(output_dir, "trajectories_all.png"), p1, 
         width = 10, height = 8, dpi = 300)
  
  # Faceted by date (sample of dates if too many)
  unique_dates <- unique(data$Date)
  
  if (length(unique_dates) > 20) {
    sample_dates <- sample(unique_dates, 20)
    data_sample <- data %>% filter(Date %in% sample_dates)
  } else {
    data_sample <- data
  }
  
  p2 <- ggplot(data_sample, aes(x = longitude, y = latitude, color = species)) +
    geom_path(aes(group = species), alpha = 0.5) +
    geom_point(size = 1) +
    scale_color_manual(values = c("Baboon" = "#D55E00", "Vervet" = "#009E73")) +
    facet_wrap(~Date, ncol = 4) +
    theme_minimal() +
    theme(axis.text = element_text(size = 6)) +
    labs(title = "Daily Movement Trajectories",
         x = "Longitude", y = "Latitude") +
    coord_fixed(ratio = 1)
  
  ggsave(file.path(output_dir, "trajectories_by_date.png"), p2, 
         width = 12, height = 12, dpi = 300)
}

# ============================================================================
# PART 2: FILTER UNREALISTIC MOVEMENTS
# ============================================================================

filter_unrealistic_movements <- function(data, max_speed_kmh = 5) {
  
  # Maximum speed: 5 km/h is reasonable for primates
  # (baboons typically 3-4 km/h, vervets similar)
  
  cat("Filtering unrealistic movements...\n")
  cat("Maximum speed threshold:", max_speed_kmh, "km/h\n\n")
  
  # Ensure data is sorted
  data <- data %>%
    arrange(species, New_Timestamp)
  
  # Calculate step lengths and speeds
  data_filtered <- data %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      # Calculate distance to next point (meters)
      next_lon = lead(longitude),
      next_lat = lead(latitude),
      next_time = lead(New_Timestamp),
      
      # Distance in meters using Haversine formula
      distance_m = distHaversine(
        cbind(longitude, latitude),
        cbind(next_lon, next_lat)
      ),
      
      # Time difference in hours
      time_diff_hours = as.numeric(difftime(next_time, New_Timestamp, units = "hours")),
      
      # Speed in km/h
      speed_kmh = (distance_m / 1000) / time_diff_hours,
      
      # Flag unrealistic movements (within same observation session)
      # Only flag if time difference is less than 2 hours (within same session)
      unrealistic = !is.na(speed_kmh) & speed_kmh > max_speed_kmh & time_diff_hours < 2
    ) %>%
    ungroup()
  
  # Summary of flagged points
  n_unrealistic <- sum(data_filtered$unrealistic, na.rm = TRUE)
  
  cat("Flagged movements:\n")
  cat(sprintf("  Total points: %d\n", nrow(data_filtered)))
  cat(sprintf("  Unrealistic movements: %d (%.2f%%)\n", 
              n_unrealistic, 100 * n_unrealistic / nrow(data_filtered)))
  
  if (n_unrealistic > 0) {
    cat("\nSpeed distribution of flagged points:\n")
    flagged <- data_filtered %>% filter(unrealistic == TRUE)
    print(summary(flagged$speed_kmh))
    
    # Plot speed distribution
    p <- ggplot(data_filtered %>% filter(!is.na(speed_kmh)), 
                aes(x = speed_kmh)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = max_speed_kmh, color = "red", 
                 linetype = "dashed", size = 1) +
      scale_x_log10() +
      theme_minimal() +
      labs(title = "Distribution of Movement Speeds",
           subtitle = paste("Red line indicates threshold:", max_speed_kmh, "km/h"),
           x = "Speed (km/h, log scale)", y = "Count")
    
    print(p)
  }
  
  # Remove unrealistic points
  data_clean <- data_filtered %>%
    filter(!unrealistic | is.na(unrealistic)) %>%
    dplyr::select(-next_lon, -next_lat, -next_time, -distance_m, 
           -time_diff_hours, -speed_kmh, -unrealistic)
  
  cat(sprintf("\nRemaining points after filtering: %d\n\n", nrow(data_clean)))
  
  return(data_clean)
}

# ============================================================================
# PART 3: COEFFICIENT OF SOCIALITY (Cs)
# ============================================================================

calculate_sociality_coefficient <- function(baboon_data, vervet_data) {
  
  cat("Calculating Coefficient of Sociality (Cs)...\n")
  cat("Using WildlifeDI package\n\n")
  
  # Prepare data for wildlifeDI
  # Need to convert to ltraj objects
  
  # Baboon trajectory
  bab_df <- baboon_data %>%
    arrange(New_Timestamp) %>%
    select(longitude, latitude, New_Timestamp) %>%
    filter(!is.na(New_Timestamp))
  
  bab_ltraj <- as.ltraj(
    xy = bab_df[, c("longitude", "latitude")],
    date = bab_df$New_Timestamp,
    id = "Baboon",
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  
  # Vervet trajectory
  ver_df <- vervet_data %>%
    arrange(New_Timestamp) %>%
    select(longitude, latitude, New_Timestamp) %>%
    filter(!is.na(New_Timestamp))
  
  ver_ltraj <- as.ltraj(
    xy = ver_df[, c("longitude", "latitude")],
    date = ver_df$New_Timestamp,
    id = "Vervet",
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  
  # Combine trajectories
  combined_ltraj <- c(bab_ltraj, ver_ltraj)
  
  # Calculate proximity events (within specified distance threshold)
  # Test multiple distance thresholds
  distance_thresholds <- c(50, 100, 150, 200, 300)  # meters
  
  proximity_results <- list()
  
  for (thresh in distance_thresholds) {
    cat(sprintf("Testing distance threshold: %d meters\n", thresh))
    
    # Calculate Cs (Coefficient of Sociality)
    # Cs ranges from -1 (avoidance) to +1 (attraction)
    prox <- GetSimultaneous(combined_ltraj, tc = 60 * 20)  # 20 minute tolerance
    
    if (!is.null(prox) && nrow(prox) > 0) {
      # Calculate distances between simultaneous fixes
      prox$distance <- distHaversine(
        cbind(prox$x1, prox$y1),
        cbind(prox$x2, prox$y2)
      )
      
      # Count proximity events
      n_close <- sum(prox$distance <= thresh, na.rm = TRUE)
      n_total <- nrow(prox)
      prop_close <- n_close / n_total
      
      proximity_results[[as.character(thresh)]] <- list(
        threshold = thresh,
        n_simultaneous = n_total,
        n_proximity = n_close,
        proportion = prop_close,
        mean_distance = mean(prox$distance, na.rm = TRUE),
        median_distance = median(prox$distance, na.rm = TRUE)
      )
      
      cat(sprintf("  Simultaneous fixes: %d\n", n_total))
      cat(sprintf("  Within %dm: %d (%.1f%%)\n", thresh, n_close, 100*prop_close))
      cat(sprintf("  Mean distance: %.1fm\n\n", mean(prox$distance, na.rm = TRUE)))
    }
  }
  
  # Calculate Coefficient of Sociality using wildlifeDI
  # This compares observed proximities to expected under random movement
  tryCatch({
    cs_result <- Cs(combined_ltraj, tc = 60 * 20, dc = 150)  # 150m threshold
    
    cat("\nCoefficient of Sociality (Cs):\n")
    cat(sprintf("  Cs value: %.3f\n", cs_result))
    cat("  Interpretation:\n")
    cat("    Cs > 0: Attraction (more proximity than expected by chance)\n")
    cat("    Cs = 0: Random association\n")
    cat("    Cs < 0: Avoidance (less proximity than expected by chance)\n\n")
    
    proximity_results$cs_value <- cs_result
  }, error = function(e) {
    cat("\nNote: Could not calculate Cs statistic.\n")
    cat("Error:", e$message, "\n\n")
  })
  
  return(list(
    proximity_results = proximity_results,
    simultaneous_fixes = prox
  ))
}

# ============================================================================
# PART 4: SSF-DIST (Step Selection Function with Distance Covariate)
# ============================================================================

analyze_ssf_dist <- function(baboon_data, vervet_data) {
  
  cat("Analyzing interactions using SSF-DIST approach...\n")
  cat("Following Muff et al. 2023 and Oliveira-Santos et al. 2022\n\n")
  
  # Prepare data for amt package
  
  # Baboon track
  bab_track <- baboon_data %>%
    arrange(New_Timestamp) %>%
    filter(!is.na(New_Timestamp)) %>%
    make_track(longitude, latitude, New_Timestamp, 
               crs = 4326, all_cols = TRUE) %>%
    transform_coords(32736)  # UTM zone 36S for Uganda
  
  # Vervet track
  ver_track <- vervet_data %>%
    arrange(New_Timestamp) %>%
    filter(!is.na(New_Timestamp)) %>%
    make_track(longitude, latitude, New_Timestamp, 
               crs = 4326, all_cols = TRUE) %>%
    transform_coords(32736)  # UTM zone 36S
  
  # Add distance to other species as covariate
  # For each baboon location, calculate distance to nearest vervet location in time
  
  # Function to get nearest location in time
  get_nearest_distance <- function(time, focal_x, focal_y, other_track) {
    time_diffs <- abs(difftime(other_track$t_, time, units = "mins"))
    nearest_idx <- which.min(time_diffs)
    
    if (length(nearest_idx) > 0 && time_diffs[nearest_idx] <= 60) {  # Within 1 hour
      dist <- sqrt((focal_x - other_track$x_[nearest_idx])^2 + 
                     (focal_y - other_track$y_[nearest_idx])^2)
      return(dist)
    } else {
      return(NA)
    }
  }
  
  # Add distance covariate to baboon track
  bab_track$dist_to_vervet <- mapply(
    get_nearest_distance,
    bab_track$t_,
    bab_track$x_,
    bab_track$y_,
    MoreArgs = list(other_track = ver_track)
  )
  
  # Add distance covariate to vervet track
  ver_track$dist_to_baboon <- mapply(
    get_nearest_distance,
    ver_track$t_,
    ver_track$x_,
    ver_track$y_,
    MoreArgs = list(other_track = bab_track)
  )
  
  # Generate random steps for SSF
  cat("Generating random steps for baboons...\n")
  
  bab_ssf_data <- bab_track %>%
    steps_by_burst() %>%
    random_steps(n_control = 10) %>%
    extract_covariates(ver_track, where = "end", name = "dist_to_vervet")
  
  cat("Generating random steps for vervets...\n")
  
  ver_ssf_data <- ver_track %>%
    steps_by_burst() %>%
    random_steps(n_control = 10) %>%
    extract_covariates(bab_track, where = "end", name = "dist_to_baboon")
  
  # Fit SSF models
  cat("\nFitting SSF models...\n\n")
  
  # Baboon SSF: Does distance to vervets affect baboon movement?
  cat("Baboon SSF (effect of vervet proximity):\n")
  
  bab_ssf_model <- tryCatch({
    fit_issf(
      case_ ~ 
        sl_ + log_sl_ +  # Step length
        cos_ta_ +         # Turn angle
        dist_to_vervet +  # Distance to vervets
        strata(step_id_),
      data = bab_ssf_data,
      model = TRUE
    )
  }, error = function(e) {
    cat("Error fitting baboon SSF:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(bab_ssf_model)) {
    print(summary(bab_ssf_model))
    cat("\nInterpretation:\n")
    cat("  Positive coefficient for dist_to_vervet: Baboons move AWAY from vervets\n")
    cat("  Negative coefficient for dist_to_vervet: Baboons move TOWARD vervets\n\n")
  }
  
  # Vervet SSF: Does distance to baboons affect vervet movement?
  cat("\nVervet SSF (effect of baboon proximity):\n")
  
  ver_ssf_model <- tryCatch({
    fit_issf(
      case_ ~ 
        sl_ + log_sl_ +    # Step length
        cos_ta_ +          # Turn angle
        dist_to_baboon +   # Distance to baboons
        strata(step_id_),
      data = ver_ssf_data,
      model = TRUE
    )
  }, error = function(e) {
    cat("Error fitting vervet SSF:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(ver_ssf_model)) {
    print(summary(ver_ssf_model))
    cat("\nInterpretation:\n")
    cat("  Positive coefficient for dist_to_baboon: Vervets move AWAY from baboons\n")
    cat("  Negative coefficient for dist_to_baboon: Vervets move TOWARD baboons\n\n")
  }
  
  # Visualize distance distributions
  dist_data <- data.frame(
    species = c(rep("Baboon", nrow(bab_track)), rep("Vervet", nrow(ver_track))),
    distance = c(bab_track$dist_to_vervet, ver_track$dist_to_baboon)
  ) %>%
    filter(!is.na(distance))
  
  p <- ggplot(dist_data, aes(x = distance, fill = species)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("Baboon" = "#D55E00", "Vervet" = "#009E73")) +
    theme_minimal() +
    labs(title = "Distribution of Inter-species Distances",
         x = "Distance to Other Species (meters)",
         y = "Count",
         fill = "Focal Species")
  
  print(p)
  
  return(list(
    baboon_model = bab_ssf_model,
    vervet_model = ver_ssf_model,
    distance_data = dist_data
  ))
}

# ============================================================================
# MAIN ANALYSIS WORKFLOW
# ============================================================================

run_full_analysis <- function(matching_data, output_dir = "movement_analysis") {
  
  cat(strrep("=", 70), "\n")
  cat("MOVEMENT ANALYSIS AND INTERACTION INFERENCE\n")
  cat(strrep("=", 70), "\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Visualize trajectories
  cat("STEP 1: Visualizing movement trajectories...\n")
  cat(strrep("-", 70), "\n")
  visualize_trajectories(matching_data, output_dir)
  
  # Step 2: Filter unrealistic movements
  cat("\nSTEP 2: Filtering unrealistic movements...\n")
  cat(strrep("-", 70), "\n")
  data_clean <- filter_unrealistic_movements(matching_data, max_speed_kmh = 5)
  
  # Separate by species
  baboon_data <- data_clean %>% filter(species == "Baboon")
  vervet_data <- data_clean %>% filter(species == "Vervet")
  
  cat(sprintf("Baboon locations: %d\n", nrow(baboon_data)))
  cat(sprintf("Vervet locations: %d\n\n", nrow(vervet_data)))
  
  # Step 3: Calculate Coefficient of Sociality
  cat("STEP 3: Calculating Coefficient of Sociality...\n")
  cat(strrep("-", 70), "\n")
  cs_results <- calculate_sociality_coefficient(baboon_data, vervet_data)
  
  # Step 4: SSF-DIST Analysis
  cat("\nSTEP 4: SSF-DIST Analysis...\n")
  cat(strrep("-", 70), "\n")
  ssf_results <- analyze_ssf_dist(baboon_data, vervet_data)
  
  # Save results
  results <- list(
    data_clean = data_clean,
    cs_results = cs_results,
    ssf_results = ssf_results
  )
  
  saveRDS(results, file.path(output_dir, "analysis_results.rds"))
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(strrep("=", 70), "\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(results)
}

# ============================================================================
# RUN ANALYSIS
# ============================================================================

# Uncomment to run full analysis:
results <- run_full_analysis(matching_data)
