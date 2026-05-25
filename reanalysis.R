
library(ggplot2)
remove.packages("rlang")
install.packages("rlang")

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

# species trajectories plots

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

# filter unrealistic movement

# Function to calculate speeds and filter (5km/hr)
filter_by_speed <- function(data, species_name, max_speed_kmh = 5) {
  
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
str(vervet_clean)
# Combine cleaned data
matching_data_clean <- bind_rows(baboon_clean, vervet_clean)

# save
saveRDS(vervet_clean, "vervet_clean.rds")
saveRDS(baboon_clean, "baboon_clean.rds")


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

#Found 5869 simultaneous fixes

#Proximity analysis:
#  Within 50m: 5 / 5869 (0.1%)
# Within 100m: 26 / 5869 (0.4%)
# Within 150m: 105 / 5869 (1.8%)
# Within 200m: 203 / 5869 (3.5%)
# Within 300m: 519 / 5869 (8.8%)

# Mean distance: 1443.4m
# Median distance: 1278.0m

# Coefficient of Sociality (simplified calculation):
#  Proportion within 150m: 0.018
# Interpretation: Low association (possible avoidance)

