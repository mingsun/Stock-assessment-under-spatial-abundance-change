library(tidyverse)
library(sf)
library(concaveman)
library(gstat)
library(terra)
library(units)

# use survey full station to define a spatial domain


# 1. summer flounder Biglow tows ----

## 1.1 Read tow data ----
station_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv") |>
  filter(DECDEG_BEGLAT <= 41.6) |>
  # filter(DECDEG_BEGLON <= -69, DECDEG_BEGLAT <= 41.6) |>
  select(YEAR, LON = DECDEG_BEGLON, LAT = DECDEG_BEGLAT) # drop rows with missing coordinates

plot(LAT ~ LON, station_df)

# Convert to sf points (WGS84)
station_sf <- st_as_sf(
  station_df,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

## ------------------------------------------------------------ ##


## 1.2 Choose a projected CRS (for 1 km grid + variogram/kriging) ----

cent <- st_coordinates(st_centroid(st_union(station_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

station_utm <- st_transform(station_sf, crs = epsg_utm)

## ------------------------------------------------------------ ##


## 1.3 Build an irregular survey-domain polygon from tow points ----

#    - concave hull, then buffer to include nearby sediment points and avoid edge artifacts

# Create a single multipoint geometry
mp <- st_union(station_utm)

pts <- st_cast(mp, "POINT", warn = FALSE)
dom <- concaveman::concaveman(pts)

# Buffer domain (meters). Adjust if you want tighter/looser.
# This helps kriging near edges and includes sediment points just outside the hull.
dom_buf <- st_buffer(dom, dist = 2 * 1000) # n km buffer

# Visual check of domain and tow points
plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(station_utm), add = TRUE, pch = 16, cex = 0.3)

# Save domain polygon for inspection
st_write(
  st_as_sf(dom_buf),
  "manuscript/3. stock assessment with spatial abundance change/data/SF_survey_domain_polygon.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

st_crs(dom_buf)$epsg

## ------------------------------------------------------------ ##

remove(list = ls())

# ---------------------------------------- #



# 2. longfin squid BTS tows ----

## 2.1 Read tow data ----
station_df <- read.csv("results/indices for assessment/squid/catch_by_tow.csv") |>
  filter(DECDEG_BEGLAT <= 41.5) |>
  # filter(DECDEG_BEGLON <= -69, DECDEG_BEGLAT <= 41.6) |>
  select(YEAR, LON = DECDEG_BEGLON, LAT = DECDEG_BEGLAT) |> # drop rows with missing coordinates
  filter(!is.na(LON))

plot(LAT ~ LON, station_df)

# Convert to sf points (WGS84)
station_sf <- st_as_sf(
  station_df,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

## ------------------------------------------------------------ ##


## 2.2 Choose a projected CRS (for 1 km grid + variogram/kriging) ----

cent <- st_coordinates(st_centroid(st_union(station_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

station_utm <- st_transform(station_sf, crs = epsg_utm)

## ------------------------------------------------------------ ##


## 2.3 Build an irregular survey-domain polygon from tow points ----

#    - concave hull, then buffer to include nearby sediment points and avoid edge artifacts

# Create a single multipoint geometry
mp <- st_union(station_utm)

pts <- st_cast(mp, "POINT", warn = FALSE)
dom <- concaveman::concaveman(pts)

# Buffer domain (meters). Adjust if you want tighter/looser.
# This helps kriging near edges and includes sediment points just outside the hull.
dom_buf <- st_buffer(dom, dist = 2 * 1000) # n km buffer

# Visual check of domain and tow points
plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(station_utm), add = TRUE, pch = 16, cex = 0.3)

# Save domain polygon for inspection
st_write(
  st_as_sf(dom_buf),
  "manuscript/3. stock assessment with spatial abundance change/data/LS_survey_domain_polygon.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

st_crs(dom_buf)$epsg

## ------------------------------------------------------------ ##

remove(list = ls())

# ---------------------------------------- #


# 3. surfclam tows: west part ----

## 3.1 Read tow data ----
station_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") |>
  filter(YEAR %in% c(2012, 2015, 2018, 2022), REGION == "SVAtoSNE") |>
  filter(LON <= -71) |>
  select(YEAR, LON, LAT) |> # drop rows with missing coordinates
  filter(!is.na(LON))

plot(LAT ~ LON, station_df)

# Convert to sf points (WGS84)
station_sf <- st_as_sf(
  station_df,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

## ------------------------------------------------------------ ##


## 3.2 Choose a projected CRS (for 1 km grid + variogram/kriging) ----

cent <- st_coordinates(st_centroid(st_union(station_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

station_utm <- st_transform(station_sf, crs = epsg_utm)

## ------------------------------------------------------------ ##


## 3.3 Build an irregular survey-domain polygon from tow points ----

#    - concave hull, then buffer to include nearby sediment points and avoid edge artifacts

# Create a single multipoint geometry
mp <- st_union(station_utm)

pts <- st_cast(mp, "POINT", warn = FALSE)
dom <- concaveman::concaveman(pts)

# Buffer domain (meters). Adjust if you want tighter/looser.
# This helps kriging near edges and includes sediment points just outside the hull.
dom_buf <- st_buffer(dom, dist = 2 * 1000) # n km buffer

# Visual check of domain and tow points
plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(station_utm), add = TRUE, pch = 16, cex = 0.3)

# Save domain polygon for inspection
st_write(
  st_as_sf(dom_buf),
  "manuscript/3. stock assessment with spatial abundance change/data/SC_west_survey_domain_polygon.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

st_crs(dom_buf)$epsg

## ------------------------------------------------------------ ##

remove(list = ls())

# ---------------------------------------- #



# 4. surfclam tows: east part ----

## 4.1 Read tow data ----
station_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") |>
  filter(YEAR %in% c(2012, 2015, 2018, 2022), REGION == "SVAtoSNE") |>
  filter(LON >= -71) |>
  select(YEAR, LON, LAT) |> # drop rows with missing coordinates
  filter(!is.na(LON))

plot(LAT ~ LON, station_df)

# Convert to sf points (WGS84)
station_sf <- st_as_sf(
  station_df,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

## ------------------------------------------------------------ ##


## 4.2 Choose a projected CRS (for 1 km grid + variogram/kriging) ----

cent <- st_coordinates(st_centroid(st_union(station_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

station_utm <- st_transform(station_sf, crs = epsg_utm)

## ------------------------------------------------------------ ##


## 4.3 Build an irregular survey-domain polygon from tow points ----

#    - concave hull, then buffer to include nearby sediment points and avoid edge artifacts

# Create a single multipoint geometry
mp <- st_union(station_utm)

pts <- st_cast(mp, "POINT", warn = FALSE)
dom <- concaveman::concaveman(pts)

# Buffer domain (meters). Adjust if you want tighter/looser.
# This helps kriging near edges and includes sediment points just outside the hull.
dom_buf <- st_buffer(dom, dist = 2 * 1000) # n km buffer

# Visual check of domain and tow points
plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(station_utm), add = TRUE, pch = 16, cex = 0.3)

# Save domain polygon for inspection
st_write(
  st_as_sf(dom_buf),
  "manuscript/3. stock assessment with spatial abundance change/data/SC_east_survey_domain_polygon.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

st_crs(dom_buf)$epsg

## ------------------------------------------------------------ ##

remove(list = ls())

# ---------------------------------------- #

