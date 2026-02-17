# Vervet cost surface
# aim is to describe areas that are too difficult for this group to climb/places they
# wont go

# Derives movement barriers from observed vervet elevation changes.
# The maximum elevation change a vervet is observed to cross in a single
# step defines the threshold beyond which terrain is impassable.


library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(amt)

# extract DEM values for each vervet location

vervet_sf <- st_as_sf(vervet_clean, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 32735)

vervet_coords <- st_coordinates(vervet_sf)

vervet_clean$elev <- extract(dem, vervet_coords)[, 1]
vervet_clean$x_utm <- vervet_coords[, 1]
vervet_clean$y_utm <- vervet_coords[, 2]

# calculate elevation change
# Sort by timestamp and calculate consecutive elevation differences
vervet_steps_elev <- vervet_clean %>%
  arrange(New_Timestamp) %>%
  mutate(
    elev_next    = lead(elev),
    time_next    = lead(New_Timestamp),
    x_next       = lead(x_utm),
    y_next       = lead(y_utm),
    time_gap_min = as.numeric(difftime(time_next, New_Timestamp, units = "mins")),
    horiz_dist   = sqrt((x_next - x_utm)^2 + (y_next - y_utm)^2),
    elev_change  = abs(elev_next - elev)
  ) %>%
  # Keep only steps within a reasonable time window (5-40 mins)
  # Excludes missing points (large gaps) that would inflate elevation change
  filter(time_gap_min >= 5, time_gap_min <= 40) %>%
  filter(!is.na(elev_change))

cat(sprintf("  %d consecutive steps within 5-40 min window\n\n", nrow(vervet_steps_elev)))

# distribution plot for deciding a threshold

# Summary statistics
elev_quantiles <- quantile(vervet_steps_elev$elev_change,
                           probs = c(0.5, 0.75, 0.9, 0.95, 0.99, 1.0))

cat("Elevation change per step (metres):\n")
cat(sprintf("  Median:      %.1f m\n", elev_quantiles["50%"]))
cat(sprintf("  75th pctile: %.1f m\n", elev_quantiles["75%"]))
cat(sprintf("  90th pctile: %.1f m\n", elev_quantiles["90%"]))
cat(sprintf("  95th pctile: %.1f m\n", elev_quantiles["95%"]))
cat(sprintf("  99th pctile: %.1f m\n", elev_quantiles["99%"]))
cat(sprintf("  Maximum:     %.1f m\n\n", elev_quantiles["100%"]))

# Full distribution plot
p1 <- ggplot(vervet_steps_elev, aes(x = elev_change)) +
  geom_histogram(bins = 60, fill = "#009E73", alpha = 0.7, color = "white") +
  geom_vline(xintercept = elev_quantiles["90%"],  linetype = "dashed", color = "blue",   size = 0.8) +
  geom_vline(xintercept = elev_quantiles["95%"],  linetype = "dashed", color = "orange", size = 0.8) +
  geom_vline(xintercept = elev_quantiles["99%"],  linetype = "dashed", color = "red",    size = 0.8) +
  annotate("text", x = elev_quantiles["90%"],  y = Inf, label = "90th", vjust = 2, hjust = -0.2, color = "blue",   size = 3.5) +
  annotate("text", x = elev_quantiles["95%"],  y = Inf, label = "95th", vjust = 2, hjust = -0.2, color = "orange", size = 3.5) +
  annotate("text", x = elev_quantiles["99%"],  y = Inf, label = "99th", vjust = 2, hjust = -0.2, color = "red",    size = 3.5) +
  theme_minimal() +
  labs(title = "Distribution of Elevation Change per Vervet Step",
       subtitle = "Dashed lines show candidate cliff barrier thresholds",
       x = "Elevation Change per Step (m)", y = "Count")

ggsave("movement_analysis/vervet_elev_change_distribution.png", p1,
       width = 9, height = 6, dpi = 300)

# Tail plot (upper 10%) to see the cliff candidates more clearly
p2 <- ggplot(vervet_steps_elev %>% filter(elev_change > elev_quantiles["75%"]),
             aes(x = elev_change)) +
  geom_histogram(bins = 40, fill = "#D55E00", alpha = 0.7, color = "white") +
  geom_vline(xintercept = elev_quantiles["90%"],  linetype = "dashed", color = "blue",   size = 0.8) +
  geom_vline(xintercept = elev_quantiles["95%"],  linetype = "dashed", color = "orange", size = 0.8) +
  geom_vline(xintercept = elev_quantiles["99%"],  linetype = "dashed", color = "red",    size = 0.8) +
  theme_minimal() +
  labs(title = "Upper Tail of Elevation Change Distribution",
       x = "Elevation Change per Step (m)", y = "Count")

ggsave("movement_analysis/vervet_elev_change_tail.png", p2,
       width = 9, height = 6, dpi = 300)

# based on viewing the upper percentiles we would say 65m and above are basically impassable

# threshold based on looking at upper bounds
CLIFF_THRESHOLD <- 30

# identify cliffs in the dem

# Calculate elevation change in every direction from each cell
# Use focal max - focal min over a 3x3 window as local relief
local_relief <- focal(dem, w = 3, fun = function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# A cell is a cliff barrier if its local relief exceeds the observed threshold
cliff_mask <- local_relief > CLIFF_THRESHOLD

n_cliff_cells <- sum(values(cliff_mask), na.rm = TRUE)
pct_cliff <- 100 * n_cliff_cells / ncell(cliff_mask)

cat(sprintf("  Cliff cells (local relief > %.1f m): %d (%.1f%% of study area)\n\n",
            CLIFF_THRESHOLD, n_cliff_cells, pct_cliff))

# cost surface using this info

# Calculate slope as base cost
slope <- terrain(dem, "slope", unit = "degrees")

# Slope cost: flat = 1, steeper = higher
slope_cost <- 1 + (slope / 45)^2

# Apply cliff barrier: very high cost where terrain exceeds vervet threshold
# Not set to Inf so the solver doesn't break, but high enough to be avoided
cost_surface <- ifel(cliff_mask, slope_cost * 1000, slope_cost)

# Normalize to 1-100
cost_min <- global(cost_surface, "min", na.rm = TRUE)[[1]]
cost_max <- global(cost_surface, "max", na.rm = TRUE)[[1]]
cost_surface_norm <- (cost_surface - cost_min) / (cost_max - cost_min) * 99 + 1

# visual

png("movement_analysis/cost_surface_empirical.png", width = 12, height = 5,
    units = "in", res = 300)
par(mfrow = c(1, 3))
plot(dem,         main = "Elevation (m)",           col = terrain.colors(100))
plot(cliff_mask,  main = sprintf("Cliff barriers\n(local relief > %.1f m)", CLIFF_THRESHOLD),
     col = c("white", "firebrick"))
plot(cost_surface_norm, main = "Cost surface (1-100)", col = terrain.colors(100))
par(mfrow = c(1, 1))
dev.off()


