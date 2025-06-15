### Load required libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(circular)
library(mgcv)
library(suncalc)
library(raster)
library(geodata)
library(gratia)
library(nasapower)
library(viridis)

### 1. Load and parse raw occurrence data
botrypus_df <- read.csv("botrypus_data.csv")

botrypus_df$eventDate <- parse_date_time(botrypus_df$eventDate, orders = c("mdy", "ymd", "dmy"))
botrypus_df <- botrypus_df %>%
  filter(!is.na(eventDate)) %>%
  filter(!year(eventDate) == 0 & !is.na(month(eventDate)) & !is.na(day(eventDate))) %>%
  mutate(eventDate = as.Date(eventDate)) %>%
  mutate(eventDate = if_else(year(eventDate) > year(Sys.Date()),
                              update(eventDate, year = year(eventDate) - 100),
                              eventDate))

### 2. Calculate day of year and filter edge cases
botrypus_df <- botrypus_df %>%
  mutate(day_of_year = as.numeric(format(eventDate, "%j"))) %>%
  filter(!(day_of_year %in% c(0, 1, 364, 365)))

### 3. Latitude binning for spatial and phenological analysis
botrypus_df <- botrypus_df %>%
  mutate(lat_bin = cut(
    latitude,
    breaks = c(-60, 10, 30, 40, 50, 60, 75),
    right = FALSE,
    include.lowest = TRUE
  ))

### 4. Compute within-bin phenological variation (SD)
botrypus_df <- botrypus_df %>%
  group_by(lat_bin) %>%
  mutate(sd_day = sd(day_of_year, na.rm = TRUE)) %>%
  ungroup()

### 5. Filter data to focus on the Americas
botrypus_df_filtered <- botrypus_df %>%
  filter(latitude >= -60 & latitude <= 90,
         longitude >= -180 & longitude <= -30,
         !is.na(sd_day),
         !is.na(lat_bin))

### 6. Plot histogram of phenology by latitude bin
ggplot(botrypus_df_filtered) +
  geom_histogram(aes(x = day_of_year, fill = lat_bin), bins = 50, position = "dodge") +
  labs(
    title = "Histograms of Botrypus virginianus Observations by Latitude Bin",
    x = "Day of Year",
    y = "Number of Observations"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  facet_wrap(~lat_bin, scales = "free_y")

### 7. Plot map of phenological standard deviation by location
world <- ne_countries(scale = "medium", returnclass = "sf")
americas <- subset(world, continent %in% c("North America", "South America"))
americas <- st_transform(americas, crs = 4326)

ggplot(data = botrypus_df_filtered) +
  geom_sf(data = americas, fill = "gray90", color = "gray30", size = 0.2) +
  geom_point(aes(x = longitude, y = latitude, color = sd_day), size = 1) +
  scale_color_viridis_c() +
  labs(
    title = "Geographic Distribution of Botrypus virginianus Observations by Date Variation",
    subtitle = "Colored by the Standard Deviation of Observation Dates",
    color = "Standard Deviation of Dates"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  coord_sf(xlim = c(-180, -30), ylim = c(-60, 90))

### 7b. Map of raw day-of-year values across the Americas
ggplot(data = botrypus_df_filtered) +
  geom_sf(data = americas, fill = "gray90", color = "gray30", size = 0.2) +
  geom_point(aes(x = longitude, y = latitude, color = day_of_year), size = 1) +
  scale_color_viridis_c(name = "Day of Year") +
  labs(
    title = "Geographic Distribution of Botrypus virginianus Observations",
    subtitle = "Colored by Day of Year",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  coord_sf(xlim = c(-180, -30), ylim = c(-60, 90))

### 8. Summary statistics by latitude bin
summary_stats <- botrypus_df_filtered %>%
  group_by(lat_bin) %>%
  summarise(
    count = n(),
    mean_day_of_year = mean(day_of_year, na.rm = TRUE),
    median_day_of_year = median(day_of_year, na.rm = TRUE),
    sd_day_of_year = sd(day_of_year, na.rm = TRUE),
    min_day_of_year = min(day_of_year, na.rm = TRUE),
    max_day_of_year = max(day_of_year, na.rm = TRUE),
    iqr_day_of_year = IQR(day_of_year, na.rm = TRUE)
  )

print(summary_stats, n = 10)
View(summary_stats)

### 9. Circular histogram by normalized month
botrypus_df_filtered <- botrypus_df_filtered %>%
  mutate(lat_bin_combined = case_when(
    latitude >= -60 & latitude < 10 ~ "-60 to 10",
    latitude >= 10 & latitude < 30 ~ "10 to 30",
    latitude >= 30 & latitude < 40 ~ "30 to 40",
    latitude >= 40 & latitude < 50 ~ "40 to 50",
    latitude >= 50 & latitude <= 70 ~ "50 to 70",
    TRUE ~ NA_character_
  ))

botrypus_df_filtered$lat_bin_combined <- factor(
  botrypus_df_filtered$lat_bin_combined,
  levels = c("-60 to 10", "10 to 30", "30 to 40", "40 to 50", "50 to 70")
)

botrypus_df_monthly <- botrypus_df_filtered %>%
  mutate(month = floor((day_of_year - 1) / 30.42) + 1)

botrypus_df_monthly_counts <- botrypus_df_monthly %>%
  group_by(lat_bin_combined, month) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(lat_bin_combined) %>%
  mutate(normalized_count = count / max(count))

botrypus_df_monthly_counts$lat_bin_combined <- factor(
  botrypus_df_monthly_counts$lat_bin_combined,
  levels = c("-60 to 10", "10 to 30", "30 to 40", "40 to 50", "50 to 70")
)

ggplot(botrypus_df_monthly_counts, aes(x = circular(month, units = "degrees", zero = 0, rotation = "counter"))) +
  geom_bar(aes(y = normalized_count), fill = "skyblue", color = "black", stat = "identity") +
  coord_polar(start = 0) +
  facet_wrap(~lat_bin_combined, ncol = 2) +
  labs(
    title = "Circular Histogram of Botrypus virginianus Occurrence by Month",
    subtitle = "Normalized by Latitude Bin",
    x = "Month", y = "Frequency"
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

### 10. GDD Calculation and Latitudinal Trend (Matches Fig. 5)

# Subset data for temperate localities (25°–65°N)
botrypus_subset <- botrypus_df_filtered %>%
  filter(latitude >= 25 & latitude <= 65) %>%
  arrange(latitude) %>%
  slice(round(seq(1, n(), length.out = 100))) %>%
  mutate(
    year = lubridate::year(eventDate),
    day_of_year = lubridate::yday(eventDate)
  )

# Define function to calculate GDD
calculate_gdd <- function(lon, lat, year, end_date, base_temp = 8) {
  start <- as.Date(paste0(year, "-01-01"))
  end <- as.Date(end_date)
  if (end < start) return(NA)
  tryCatch({
    temp_series <- get_power(
      community = "ag",
      lonlat = c(lon, lat),
      pars = "T2M",
      dates = c(as.character(start), as.character(end)),
      temporal_api = "daily"
    )
    gdd <- sum(pmax(0, temp_series$T2M - base_temp), na.rm = TRUE)
    return(gdd)
  }, error = function(e) {
    message("Failed for: ", lon, ", ", lat, ", ", year)
    return(NA)
  })
}

# Calculate GDD for subset ## this takes a while and will fail for some pts.
botrypus_subset$gdd_8 <- mapply(
  calculate_gdd,
  lon = botrypus_subset$longitude,
  lat = botrypus_subset$latitude,
  year = botrypus_subset$year,
  end_date = botrypus_subset$eventDate
)

# Plot GDD vs Day of Year
ggplot(botrypus_subset, aes(x = day_of_year, y = gdd_8, color = latitude)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(name = "Latitude (°N)", option = "D", direction = 1) +
  labs(
    x = "Day of Year",
    y = "GDD to Emergence (base 8 °C)",
    title = "GDD vs. Day of Year (colored by latitude)"
  ) +
  theme_minimal()

# Plot GDD vs Latitude with LOESS fit
ggplot(botrypus_subset, aes(x = latitude, y = gdd_8)) +
  geom_point(size = 3, alpha = 0.8, color = "black", fill = "lightblue", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Latitude (°N)",
    y = "GDD to Emergence (base 8 °C)",
    title = "GDD to Emergence vs. Latitude (25°–65°N)"
  ) +
  theme_minimal()

# Linear model for GDD ~ Latitude
gdd_lm <- lm(gdd_8 ~ latitude, data = botrypus_subset)
summary(gdd_lm)
