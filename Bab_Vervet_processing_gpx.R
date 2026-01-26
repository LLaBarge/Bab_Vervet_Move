# GPX Files - Baboon and Vervet Movement Data
 
# Install packages
required_packages <- c("sf", "dplyr", "lubridate")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


gpx_folder <- "C:/Users/lrlab/PPP scan points/Baboon_Vervet_Move/gpx"

# ============================================================================
# function to read a gpx file with human-written notes
# ============================================================================

read_gpx_file <- function(gpx_path) {
  
  # Extract species from filename (Baboon or Vervet)
  filename <- basename(gpx_path)
  species <- ifelse(grepl("Baboon", filename, ignore.case = TRUE), "Baboon",
                    ifelse(grepl("Vervet", filename, ignore.case = TRUE), "Vervet", NA))
  
  if (is.na(species)) {
    warning(sprintf("Could not determine species from filename: %s", filename))
    return(NULL)
  }
  
  tryCatch({
    # Read waypoints
    layers <- st_layers(gpx_path)
    
    if ("waypoints" %in% layers$name) {
      waypoints <- st_read(gpx_path, layer = "waypoints", quiet = TRUE)
      
      if (nrow(waypoints) > 0) {
        # Extract coordinates
        coords <- st_coordinates(waypoints)
        
        # Create dataframe
        df <- data.frame(
          name = character(nrow(waypoints)),
          longitude = coords[, "X"],
          latitude = coords[, "Y"],
          species = species,
          Timestamp = character(nrow(waypoints)),
          stringsAsFactors = FALSE
        )
        
        # Extract waypoint name
        name_cols <- c("name", "Name", "NAME", "wpt_name", "desc", "cmt")
        for (col in name_cols) {
          if (col %in% names(waypoints)) {
            df$name <- as.character(waypoints[[col]])
            break
          }
        }
        
        # If no name column found, use default
        if (all(df$name == "")) {
          df$name <- paste0("waypoint_", 1:nrow(df))
        }
        
        # Extract timestamp
        time_cols <- c("time", "Time", "TIME", "timestamp", "date")
        for (col in time_cols) {
          if (col %in% names(waypoints)) {
            df$Timestamp <- as.POSIXct(waypoints[[col]])
            df$Timestamp <- format(df$Timestamp, "%Y-%m-%d %H:%M:%S")
            break
          }
        }
        
        # If no time found, set to NA
        if (all(df$Timestamp == "")) {
          df$Timestamp <- NA
        }
        
        return(df)
        
      } else {
        return(NULL)
      }
      
    } else {
      warning(sprintf("No waypoints layer found in %s", filename))
      return(NULL)
    }
    
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", filename, e$message))
    return(NULL)
  })
}

# Get all GPX files
gpx_files <- list.files(gpx_folder, pattern = "\\.gpx$", 
                        full.names = TRUE, ignore.case = TRUE)

if (length(gpx_files) == 0) {
  stop("No GPX files found in the specified folder")
}

# Read all GPX files
all_data <- list()

for (gpx_file in gpx_files) {
  result <- read_gpx_file(gpx_file)
  
  if (!is.null(result) && nrow(result) > 0) {
    all_data[[length(all_data) + 1]] <- result
  }
}

if (length(all_data) == 0) {
  stop("No data was successfully extracted from any GPX files!")
}

# Combine all data
combined_data <- bind_rows(all_data)

# ============================================================================
# SPLIT INTO AD_LIB AND WAYPOINTS
# ============================================================================

# Function to check if name starts with a digit (waypoints) or letter (ad_lib)
starts_with_digit <- function(name) {
  grepl("^\\d", trimws(name))
}

# Split data: waypoints start with digit, ad_lib starts with letter
waypoints <- combined_data %>%
  filter(starts_with_digit(name))

ad_lib <- combined_data %>%
  filter(!starts_with_digit(name))

# Create New_Timestamp column for waypoints
if (nrow(waypoints) > 0) {
  
  waypoints$New_Timestamp <- sapply(waypoints$name, function(name_str) {
    tryCatch({
      # Clean up the name: remove trailing V/B, extra spaces
      cleaned <- trimws(toupper(name_str))
      cleaned <- sub("[VB]$", "", cleaned)  # Remove trailing V or B
      cleaned <- trimws(cleaned)
      
      # Remove all spaces to make parsing easier
      no_spaces <- gsub("\\s+", "", cleaned)
      
      # Extract time (first 4 digits)
      time_match <- regexpr("^\\d{4}", no_spaces)
      if (time_match == -1) return(NA)
      
      time_str <- substr(no_spaces, 1, 4)
      hour <- substr(time_str, 1, 2)
      minute <- substr(time_str, 3, 4)
      time_formatted <- sprintf("%s:%s:00", hour, minute)
      
      # Remove the time part
      remaining <- substr(no_spaces, 5, nchar(no_spaces))
      
      # Try to extract date pattern: 2 digits + 3-4 letters + 2 digits
      # Handles both "09JUN16" and "07JUNE19"
      date_pattern <- regexpr("\\d{2}[A-Z]{3,4}\\d{2}", remaining)
      if (date_pattern == -1) return(NA)
      
      date_str <- substr(remaining, date_pattern, 
                         date_pattern + attr(date_pattern, "match.length") - 1)
      
      # Parse date components
      day <- substr(date_str, 1, 2)
      
      # Month could be 3 or 4 letters (JUN or JUNE)
      # Year is always the last 2 digits
      year <- substr(date_str, nchar(date_str) - 1, nchar(date_str))
      month_str <- substr(date_str, 3, nchar(date_str) - 2)
      
      # Convert month to number (handles both JUN and JUNE)
      # Try 3-letter match first
      month_num <- match(substr(month_str, 1, 3), toupper(month.abb))
      
      if (is.na(month_num)) return(NA)
      
      # All years are 2000s (2014-2018 range)
      full_year <- paste0("20", year)
      
      # Combine into yyyy-mm-dd hh:mm:ss format
      datetime_str <- sprintf("%s-%02d-%s %s", 
                              full_year, month_num, day, time_formatted)
      
      return(datetime_str)
      
    }, error = function(e) {
      return(NA)
    })
  })
}

# Fix specific known errors
if (nrow(waypoints) >= 1573) {
  waypoints$New_Timestamp[1573] <- sub("^2029", "2019", waypoints$New_Timestamp[1573])
}

if (nrow(waypoints) >= 13621) {
  waypoints$New_Timestamp[13621] <- sub("^2001", "2014", waypoints$New_Timestamp[13621])
}

if (nrow(waypoints) >= 484) {
  waypoints$New_Timestamp[484] <- "2015-04-11 06:20:00"
}

if (nrow(waypoints) >= 7735) {
  waypoints$New_Timestamp[7735] <- "2015-12-11 13:40:00"
}


View(waypoints)


# Create Baboon dataframe from waypoints only
babs_all <- waypoints %>%
  filter(species == "Baboon") %>%
  select(name, longitude, latitude, species, New_Timestamp)

# Create Vervet dataframe from waypoints only
ver_all <- waypoints %>%
  filter(species == "Vervet") %>%
  select(name, longitude, latitude, species, New_Timestamp)


# Remove duplicates based on exact Timestamp match
babs <- babs_all %>%
  distinct(New_Timestamp, .keep_all = TRUE) %>%
  arrange(New_Timestamp)
View(babs)

ver <- ver_all %>%
  distinct(New_Timestamp, .keep_all = TRUE) %>%
  arrange(New_Timestamp)
View(ver)

# remove na
ver_edited <- na.omit(ver)
bab_edited <- na.omit(babs)
str(ver_edited)
str(ver)
str(bab_edited)
str(babs)

# Extract date from Timestamp
bab_edited$Date <- substr(bab_edited$New_Timestamp, 1, 10)
ver_edited$Date <- substr(ver_edited$New_Timestamp, 1, 10)

# Find dates that exist in both datasets
matching_dates <- intersect(unique(bab_edited$Date), unique(ver_edited$Date))

# Create combined dataframe with rows from matching dates
matching_data <- bind_rows(
  bab_edited %>% filter(Date %in% matching_dates),
  ver_edited %>% filter(Date %in% matching_dates)
) %>%
  arrange(Date, New_Timestamp)

View(matching_data)



