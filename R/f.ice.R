# function to extract sea ice raster data to xy coordinates at time point t
# sea ice input is a raster stack with 20 layers - 4 season * 5 years
# function input is season (1:4) and year (1:5) from which time point t is calculated (1:20)

f.ice <- function(x, y, season, year) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(ice))
  proj4string(spp) <- fm_sp_get_crs(ice)
  
  # convert season and year to indices for raster stack 
  if (year == 1){
    t <- season
  } else
  if (year == 2){
    t <- 4 + season
  } else
    if (year == 3){
      t <- 8 + season
    } else
  if (year == 4){
    t <- 12 + season
  } else
  if (year == 5){
    t <- 16 + season
  }
  
  # Subset raster stack at correct time index
  ice_sub <- raster::subset(ice, t)
  # Convert to SpatialPixelsDataFrame for inlabru extract
  ice_sub <- as(ice_sub, "SpatialPixelsDataFrame")
  
  # Extract season and year specific sea ice values at spp coords
  v <- over(spp, ice_sub)
  names(v) <- "ice"
  if (any(is.na(v$ice))) {
    v$ice <- inlabru:::bru_fill_missing(ice_sub, spp, v$ice)
  }
  return(v$ice)
}

# ends
