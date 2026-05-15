# proximity and detection analysis
# Examines proximity events and behavioral response

library(sf)
library(dplyr)
library(ggplot2)
library(lubridate)

MAX_PROXIMITY_M <- 350
SCAN_SAMPLE_WINDOW_HOURS <- 1
# data load

vervet_clean <- readRDS("vervet_clean.rds")
baboon_clean <- readRDS("baboon_clean.rds")

cat("  Vervet fixes:", nrow(vervet_clean), "\n")
cat("  Baboon fixes:", nrow(baboon_clean), "\n\n")

# get column names
find_col <- function(df, patterns) {
  for (p in patterns) {
    matches <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  return(NA)
}

ver_time <- find_col(vervet_clean, c("timestamp", "New_Timestamp", "time", "t_"))
ver_lon <- find_col(vervet_clean, c("longitude", "lon"))
ver_lat <- find_col(vervet_clean, c("latitude", "lat"))

bab_time <- find_col(baboon_clean, c("timestamp", "New_Timestamp", "time", "t_"))
bab_lon <- find_col(baboon_clean, c("longitude", "lon"))
bab_lat <- find_col(baboon_clean, c("latitude", "lat"))

# ============================================================================
# ADD UTM COORDINATES 

if (!("x_utm" %in% names(vervet_clean))) {
  
  
  vervet_sf <- st_as_sf(vervet_clean, coords = c(ver_lon, ver_lat), crs = 4326) %>%
    st_transform(32735)
  
  baboon_sf <- st_as_sf(baboon_clean, coords = c(bab_lon, bab_lat), crs = 4326) %>%
    st_transform(32735)
  
  # Extract UTM coordinates
  ver_coords <- st_coordinates(vervet_sf)
  bab_coords <- st_coordinates(baboon_sf)
  
  vervet_clean$x_utm <- ver_coords[, 1]
  vervet_clean$y_utm <- ver_coords[, 2]
  
  baboon_clean$x_utm <- bab_coords[, 1]
  baboon_clean$y_utm <- bab_coords[, 2]
}

vervet_clean$timestamp_clean <- vervet_clean[[ver_time]]
baboon_clean$timestamp_clean <- baboon_clean[[bab_time]]

# Get GPS dates for reference
gps_dates <- unique(as.Date(vervet_clean$timestamp_clean, tz = "UTC"))
cat("GPS dates available:", length(gps_dates), "\n")
cat("  Range:", as.character(range(gps_dates)), "\n\n")

# ============================================================================
# find proximity events

proximity <- data.frame()

n_ver <- nrow(vervet_clean)
n_bab <- nrow(baboon_clean)

for (i in 1:n_ver) {
  
  ver_time <- vervet_clean$timestamp_clean[i]
  ver_x <- vervet_clean$x_utm[i]
  ver_y <- vervet_clean$y_utm[i]
  
  # Find baboon locations within 30 minutes
  time_diffs <- abs(difftime(baboon_clean$timestamp_clean, ver_time, units = "mins"))
  
  # Check if we have any matches
  if (length(time_diffs) == 0) next
  
  within_30min <- which(time_diffs <= 30)
  
  if (length(within_30min) == 0) next
  
  # For each nearby baboon location, calculate distance
  for (j in within_30min) {
    
    bab_x <- baboon_clean$x_utm[j]
    bab_y <- baboon_clean$y_utm[j]
    
    dist <- sqrt((ver_x - bab_x)^2 + (ver_y - bab_y)^2)
    
    if (dist <= MAX_PROXIMITY_M) {
      proximity <- rbind(proximity, data.frame(
        vervet_idx = i,
        baboon_idx = j,
        timestamp = ver_time,
        distance_m = dist,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (i %% 500 == 0) {
    cat("  Processed", i, "/", n_ver, "vervet fixes (", 
        nrow(proximity), "proximity obs so far)\r")
  }
}

cat("\nFound", nrow(proximity), "proximity observations\n\n")


# ============================================================================
# GROUP INTO ENCOUNTER EVENTS

proximity <- proximity %>% arrange(timestamp)

proximity$encounter_id <- 1

if (nrow(proximity) > 1) {
  for (i in 2:nrow(proximity)) {
    gap <- difftime(proximity$timestamp[i], proximity$timestamp[i-1], units = "mins")
    if (gap > 60) {
      proximity$encounter_id[i] <- proximity$encounter_id[i-1] + 1
    } else {
      proximity$encounter_id[i] <- proximity$encounter_id[i-1]
    }
  }
}

n_encounters <- length(unique(proximity$encounter_id))
cat("  Unique encounters:", n_encounters, "\n\n")

# ============================================================================
# ENCOUNTER SUMMARY STATISTICS

encounter_summary <- proximity %>%
  group_by(encounter_id) %>%
  summarize(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_min = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins")),
    n_obs = n(),
    min_distance_m = min(distance_m),
    mean_distance_m = mean(distance_m),
    max_distance_m = max(distance_m),
    .groups = "drop"
  )


cat("  Minimum:", round(min(encounter_summary$min_distance_m)), "m\n")
cat("  Q1:", round(quantile(encounter_summary$min_distance_m, 0.25)), "m\n")
cat("  Median:", round(median(encounter_summary$min_distance_m)), "m\n")
cat("  Q3:", round(quantile(encounter_summary$min_distance_m, 0.75)), "m\n")
cat("  Maximum:", round(max(encounter_summary$min_distance_m)), "m\n")
cat("  Mean:", round(mean(encounter_summary$min_distance_m)), "m\n\n")

# ============================================================================
# MOVEMENT DIRECTION ANALYSIS
movement <- data.frame()

for (enc_id in unique(proximity$encounter_id)) {
  
  enc_data <- proximity %>% filter(encounter_id == enc_id) %>% arrange(timestamp)
  
  if (nrow(enc_data) < 2) next
  
  for (i in 1:(nrow(enc_data)-1)) {
    
    ver_i1 <- enc_data$vervet_idx[i]
    ver_i2 <- enc_data$vervet_idx[i+1]
    
    # Skip if same vervet position
    if (ver_i1 == ver_i2) next
    
    dist1 <- enc_data$distance_m[i]
    dist2 <- enc_data$distance_m[i+1]
    
    # Check if vervet moved
    ver_x1 <- vervet_clean$x_utm[ver_i1]
    ver_y1 <- vervet_clean$y_utm[ver_i1]
    ver_x2 <- vervet_clean$x_utm[ver_i2]
    ver_y2 <- vervet_clean$y_utm[ver_i2]
    
    move_dist <- sqrt((ver_x2 - ver_x1)^2 + (ver_y2 - ver_y1)^2)
    
    # Only count if vervet actually moved
    if (move_dist > 5) {
      
      moved_away <- dist2 > dist1
      
      movement <- rbind(movement, data.frame(
        encounter_id = enc_id,
        distance_to_baboon_m = dist1,
        moved_away = moved_away,
        distance_change_m = dist2 - dist1,
        movement_dist_m = move_dist,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# ============================================================================
# detection distance
# ============================================================================

if (nrow(movement) == 0) {
  cat("No movement data available\n")} else {
  
  bins <- seq(0, MAX_PROXIMITY_M, by = 100)
  
  fid <- data.frame()
  
  for (i in 1:(length(bins)-1)) {
    
    in_bin <- movement %>%
      filter(distance_to_baboon_m >= bins[i] & distance_to_baboon_m < bins[i+1])
    
    if (nrow(in_bin) > 0) {
      
      n_away <- sum(in_bin$moved_away)
      pct_away <- 100 * n_away / nrow(in_bin)
      
      fid <- rbind(fid, data.frame(
        distance_bin_mid = (bins[i] + bins[i+1]) / 2,
        n_movements = nrow(in_bin),
        n_away = n_away,
        pct_away = pct_away
      ))
    }
  }}
  
  
  print(fid)

  # Estimate detection distance (maximum distance where >50% movements are away)
  fid_50 <- fid %>% filter(pct_away > 50)
  
  if (nrow(fid_50) > 0) {
    detection_dist <- max(fid_50$distance_bin_mid)
    cat("Detection distance:", round(detection_dist), "m\n")
    cat("(Maximum distance where >50% of movements are away from baboons)\n")
  }
}


# ============================================================================
# GPS-INFORMED DATE PARSER
# ============================================================================

parse_date_gps_informed <- function(date_str, gps_dates_ref) {
  date_str <- trimws(as.character(date_str))
  
  # NUMERIC DATES (R format)
  if (grepl("^\\d+$", date_str)) {
    r_numeric <- as.numeric(date_str)
    date <- as.Date(r_numeric, origin = "1970-01-01")
    year <- as.integer(format(date, "%Y"))
    if (!is.na(date) && year >= 2014 && year <= 2019) {
      return(date)
    }
    return(NA)
  }
  
  # TEXT DATES
  if (grepl("/", date_str)) {
    parts <- strsplit(date_str, "/")[[1]]
    if (length(parts) != 3) return(NA)
    
    part1 <- as.integer(parts[1])
    part2 <- as.integer(parts[2])
    year <- as.integer(parts[3])
    
    # UNAMBIGUOUS CASES
    if (part1 > 12) {
      return(as.Date(paste(year, part2, part1, sep = "-")))
    }
    
    if (part2 > 12) {
      return(as.Date(paste(year, part1, part2, sep = "-")))
    }
    
    # AMBIGUOUS CASE
    uk_date <- tryCatch({
      as.Date(paste(year, part2, part1, sep = "-"))
    }, error = function(e) NA)
    
    us_date <- tryCatch({
      as.Date(paste(year, part1, part2, sep = "-"))
    }, error = function(e) NA)
    
    # Check which one is in GPS data
    uk_in_gps <- !is.na(uk_date) && (uk_date %in% gps_dates_ref)
    us_in_gps <- !is.na(us_date) && (us_date %in% gps_dates_ref)
    
    if (uk_in_gps && !us_in_gps) {
      return(uk_date)
    } else if (us_in_gps && !uk_in_gps) {
      return(us_date)
    } else if (uk_in_gps && us_in_gps) {
      return(uk_date)
    } else {
      if (!is.na(uk_date)) {
        return(uk_date)
      } else {
        return(us_date)
      }
    }
  }
  
  return(NA)
}

# ============================================================================
# load ad lib data

adlib <- read.csv("Vervet_adlib.csv", stringsAsFactors = FALSE)

obs_col <- names(adlib)[grepl("OBSERVATION", names(adlib), ignore.case = TRUE)][1]
date_col <- names(adlib)[grepl("DATE", names(adlib), ignore.case = TRUE)][1]

# Find baboon observations
baboon_obs <- adlib[grepl("baboon", adlib[[obs_col]], ignore.case = TRUE), ]

# GPS-informed date parsing
baboon_obs$date_clean <- sapply(baboon_obs[[date_col]], 
                                parse_date_gps_informed, 
                                gps_dates_ref = gps_dates)
baboon_obs$date_clean <- as.Date(baboon_obs$date_clean, origin = "1970-01-01")

adlib_dates <- unique(baboon_obs$date_clean[!is.na(baboon_obs$date_clean)])

cat("  Adlib baboon observations:", nrow(baboon_obs), "\n")
cat("  Parsed dates:", sum(!is.na(baboon_obs$date_clean)), "/", nrow(baboon_obs), "\n")
cat("  Unique dates:", length(adlib_dates), "\n\n")

# Label encounters
proximity$date <- as.Date(proximity$timestamp)
proximity$explicit <- proximity$date %in% adlib_dates

cat("Encounter labels:\n")
cat("  Explicit (observed):", sum(proximity$explicit), "observations\n")
cat("  Implicit (not observed):", sum(!proximity$explicit), "observations\n\n")

# ============================================================================
# LOAD SCAN DATA
# ============================================================================

scans <- read.csv("Scans_combined.csv", stringsAsFactors = FALSE, quote = "\"")

scan_date_col <- names(scans)[grepl("Date", names(scans), ignore.case = TRUE)][1]
scan_time_col <- names(scans)[grepl("Time", names(scans), ignore.case = TRUE)][1]
scan_height_col <- names(scans)[grepl("Height", names(scans), ignore.case = TRUE)][1]

# Clean height data
clean_height <- function(h) {
  h <- trimws(as.character(h))
  
  if (toupper(h) == "UN" || h == "" || is.na(h)) {
    return(NA)
  }
  
  if (grepl("^>\\s*10", h, ignore.case = TRUE)) {
    return(11)
  }
  
  h_clean <- gsub("[^0-9.]", "", h)
  height <- as.numeric(h_clean)
  
  if (!is.na(height) && height > 10) {
    return(11)
  }
  
  return(height)
}

scans$height_clean <- sapply(scans[[scan_height_col]], clean_height)

# Parse dates
scans$date_clean <- sapply(scans[[scan_date_col]], 
                           parse_date_gps_informed, 
                           gps_dates_ref = gps_dates)
scans$date_clean <- as.Date(scans$date_clean, origin = "1970-01-01")

# Parse time
parse_scan_time <- function(time_str) {
  time_str <- trimws(as.character(time_str))
  
  if (grepl("^\\d{1,2}:\\d{2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:%02d", 
                   as.integer(parts[1]), 
                   as.integer(parts[2]), 
                   as.integer(parts[3])))
  }
  
  if (grepl("^\\d{1,2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:00", 
                   as.integer(parts[1]), 
                   as.integer(parts[2])))
  }
  
  return(NA)
}

scans$time_clean <- sapply(scans[[scan_time_col]], parse_scan_time)
scans$datetime <- as.POSIXct(paste(scans$date_clean, scans$time_clean), 
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

cat("  Scan observations:", nrow(scans), "\n")
cat("  Valid datetimes:", sum(!is.na(scans$datetime)), "\n\n")

# Create scan samples
scan_samples <- scans %>%
  filter(!is.na(datetime)) %>%
  group_by(date_clean, time_clean, datetime) %>%
  summarize(
    n_individuals = n(),
    n_valid_heights = sum(!is.na(height_clean)),
    mean_height_m = ifelse(n_valid_heights > 0, 
                           mean(height_clean, na.rm = TRUE), 
                           NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_height_m)) %>%
  arrange(datetime)

# encounter analysis

encounter_ids <- unique(proximity$encounter_id)
encounter_results <- data.frame()

for (enc_id in encounter_ids) {
  
  enc_data <- proximity %>% filter(encounter_id == enc_id) %>% arrange(timestamp)
  
  if (nrow(enc_data) < 2) next
  
  # Basic info
  start_time <- min(enc_data$timestamp)
  end_time <- max(enc_data$timestamp)
  duration_min <- as.numeric(difftime(end_time, start_time, units = "mins"))
  min_dist <- min(enc_data$distance_m)
  mean_dist <- mean(enc_data$distance_m)
  date <- as.Date(start_time)
  is_explicit <- any(enc_data$explicit)
  
  # Observation text
  obs_text <- ""
  if (is_explicit) {
    obs_match <- baboon_obs %>% filter(date_clean == date)
    if (nrow(obs_match) > 0) {
      obs_text <- paste(obs_match[[obs_col]], collapse = "; ")
    }
  }
  
  # Height change
  height_before <- NA
  height_after <- NA
  height_change <- NA
  time_before <- NA
  time_after <- NA
  
  scans_before <- scan_samples %>%
    filter(datetime < start_time,
           datetime >= (start_time - hours(SCAN_SAMPLE_WINDOW_HOURS)))
  
  if (nrow(scans_before) > 0) {
    closest_idx <- which.min(abs(difftime(scans_before$datetime, start_time, units = "mins")))
    height_before <- scans_before$mean_height_m[closest_idx]
    time_before <- scans_before$datetime[closest_idx]
  }
  
  scans_after <- scan_samples %>%
    filter(datetime > end_time,
           datetime <= (end_time + hours(SCAN_SAMPLE_WINDOW_HOURS)))
  
  if (nrow(scans_after) > 0) {
    closest_idx <- which.min(abs(difftime(scans_after$datetime, end_time, units = "mins")))
    height_after <- scans_after$mean_height_m[closest_idx]
    time_after <- scans_after$datetime[closest_idx]
  }
  
  if (!is.na(height_before) && !is.na(height_after)) {
    height_change <- height_after - height_before
  }
  
  # Movement analysis
  movements_toward <- 0
  movements_away <- 0
  step_length_toward <- c()
  step_length_away <- c()
  angular_deviations <- c()
  
  for (i in 1:(nrow(enc_data)-1)) {
    
    ver_i1 <- enc_data$vervet_idx[i]
    ver_i2 <- enc_data$vervet_idx[i+1]
    
    if (ver_i1 == ver_i2) next
    
    ver_x1 <- vervet_clean$x_utm[ver_i1]
    ver_y1 <- vervet_clean$y_utm[ver_i1]
    ver_x2 <- vervet_clean$x_utm[ver_i2]
    ver_y2 <- vervet_clean$y_utm[ver_i2]
    
    bab_i <- enc_data$baboon_idx[i]
    bab_x <- baboon_clean$x_utm[bab_i]
    bab_y <- baboon_clean$y_utm[bab_i]
    
    dist_before <- sqrt((ver_x1 - bab_x)^2 + (ver_y1 - bab_y)^2)
    dist_after <- sqrt((ver_x2 - bab_x)^2 + (ver_y2 - bab_y)^2)
    
    ver_dx <- ver_x2 - ver_x1
    ver_dy <- ver_y2 - ver_y1
    move_dist <- sqrt(ver_dx^2 + ver_dy^2)
    
    if (move_dist < 10) next
    
    to_bab_x <- bab_x - ver_x1
    to_bab_y <- bab_y - ver_y1
    to_bab_dist <- sqrt(to_bab_x^2 + to_bab_y^2)
    
    dot_product <- ver_dx * to_bab_x + ver_dy * to_bab_y
    cos_angle <- dot_product / (move_dist * to_bab_dist)
    cos_angle <- max(-1, min(1, cos_angle))
    angle_rad <- acos(cos_angle)
    angle_deg <- angle_rad * 180 / pi
    
    angular_deviation <- angle_deg
    
    if (dist_after > dist_before) {
      movements_away <- movements_away + 1
      step_length_away <- c(step_length_away, move_dist)
    } else {
      movements_toward <- movements_toward + 1
      step_length_toward <- c(step_length_toward, move_dist)
    }
    
    angular_deviations <- c(angular_deviations, angular_deviation)
  }
  
  # Summary stats
  total_movements <- movements_toward + movements_away
  pct_away <- if (total_movements > 0) 100 * movements_away / total_movements else NA
  
  mean_step_toward <- if (length(step_length_toward) > 0) mean(step_length_toward) else NA
  mean_step_away <- if (length(step_length_away) > 0) mean(step_length_away) else NA
  step_length_diff <- if (!is.na(mean_step_away) && !is.na(mean_step_toward)) {
    mean_step_away - mean_step_toward
  } else NA
  
  mean_angle <- if (length(angular_deviations) > 0) mean(angular_deviations) else NA
  
  # Store results
  encounter_results <- rbind(encounter_results, data.frame(
    encounter_id = enc_id,
    start_time = start_time,
    end_time = end_time,
    duration_min = duration_min,
    date = date,
    n_observations = nrow(enc_data),
    min_distance_m = min_dist,
    mean_distance_m = mean_dist,
    detection_type = if (is_explicit) "explicit" else "implicit",
    observation = obs_text,
    height_before_m = height_before,
    height_after_m = height_after,
    height_change_m = height_change,
    scan_time_before = time_before,
    scan_time_after = time_after,
    total_movements = total_movements,
    movements_toward = movements_toward,
    movements_away = movements_away,
    pct_movements_away = pct_away,
    mean_step_toward_m = mean_step_toward,
    mean_step_away_m = mean_step_away,
    step_length_difference_m = step_length_diff,
    mean_angular_deviation_deg = mean_angle,
    stringsAsFactors = FALSE
  ))
  
  if (enc_id %% 10 == 0) cat("  Processed", enc_id, "/", length(encounter_ids), "encounters\r")
}



# ============================================================================
# Summary stats
cat("Total encounters:", nrow(encounter_results), "\n")
cat("  Explicit (observed):", sum(encounter_results$detection_type == "explicit"), "\n")
cat("  Implicit (not observed):", sum(encounter_results$detection_type == "implicit"), "\n\n")

cat("  Minimum:", round(min(encounter_results$min_distance_m)), "m\n")
cat("  Median:", round(median(encounter_results$min_distance_m)), "m\n")
cat("  Mean:", round(mean(encounter_results$min_distance_m)), "m\n\n")

# Height change
n_height <- sum(!is.na(encounter_results$height_change_m))
if (n_height > 0) {
  cat("Height change from scan samples (n =", n_height, "):\n")
  cat("  Mean:", sprintf("%.2f m", mean(encounter_results$height_change_m, na.rm=TRUE)), "\n")
  cat("  Median:", sprintf("%.2f m", median(encounter_results$height_change_m, na.rm=TRUE)), "\n")
  cat("  Increased:", sum(encounter_results$height_change_m > 0, na.rm=TRUE), "\n")
  cat("  Decreased:", sum(encounter_results$height_change_m < 0, na.rm=TRUE), "\n\n")
}

# Movement
n_move <- sum(!is.na(encounter_results$pct_movements_away))
if (n_move > 0) {
  cat("Movement away (n =", n_move, "encounters):\n")
  cat("  Mean % away:", sprintf("%.1f%%", mean(encounter_results$pct_movements_away, na.rm=TRUE)), "\n")
  cat("  Encounters mostly away (>50%):", 
      sum(encounter_results$pct_movements_away > 50, na.rm=TRUE), "\n\n")
}

# Explicit vs implicit
if (nrow(encounter_results) > 0) {
  explicit_data <- encounter_results %>% filter(detection_type == "explicit")
  implicit_data <- encounter_results %>% filter(detection_type == "implicit")
  
  if (nrow(explicit_data) > 0 && nrow(implicit_data) > 0) {
    cat("Explicit vs Implicit encounters:\n")
    cat("  Minimum distance:\n")
    cat("    Explicit:", sprintf("%.0f m", mean(explicit_data$min_distance_m, na.rm=TRUE)), "\n")
    cat("    Implicit:", sprintf("%.0f m", mean(implicit_data$min_distance_m, na.rm=TRUE)), "\n")
    cat("  % movements away:\n")
    cat("    Explicit:", sprintf("%.1f%%", mean(explicit_data$pct_movements_away, na.rm=TRUE)), "\n")
    cat("    Implicit:", sprintf("%.1f%%", mean(implicit_data$pct_movements_away, na.rm=TRUE)), "\n\n")
  }
}


# CSV files
write.csv(proximity, "all_proximity_observations.csv", row.names = FALSE)
write.csv(encounter_summary, "encounter_summary.csv", row.names = FALSE)
write.csv(encounter_results, "encounter_behavioral_analysis.csv", row.names = FALSE)

if (nrow(movement) > 0) {
  write.csv(movement, "movement_analysis.csv", row.names = FALSE)
}

if (exists("fid") && nrow(fid) > 0) {
  write.csv(fid, "detection_events.csv", row.names = FALSE)
}

# RDS files for further analysis
saveRDS(proximity, "proximity.rds")
saveRDS(encounter_summary, "encounter_summary.rds")
saveRDS(encounter_results, "encounter_results.rds")

# plots
p1 <- ggplot(encounter_summary, aes(x = min_distance_m)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Closest Approaches",
       x = "Minimum Distance (m)", y = "Count") +
  geom_vline(xintercept = median(encounter_summary$min_distance_m),
             color = "red", linetype = "dashed", linewidth = 1)

ggsave("closest_approach_distribution.png", p1, width = 8, height = 5, dpi = 300)


if (exists("fid") && nrow(fid) > 0) {
  p2 <- ggplot(fid, aes(x = distance_bin_mid, y = pct_away)) +
    geom_line(color = "darkgreen", linewidth = 1) +
    geom_point(aes(size = n_movements), color = "darkgreen") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
         x = "Distance to Baboons (m)",
         y = "% Movements Away",
         size = "N movements")
  
  ggsave("flight_initiation_distance.png", p2, width = 8, height = 5, dpi = 300)
}
p2

