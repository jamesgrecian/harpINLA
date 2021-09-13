# function to extract bathymetry data to xy coordinates

f.bat <- function(x, y) {
  
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(coords = cbind(x = x, y = y), proj4string = fm_sp_get_crs(bat_proj))
    
  # Convert to SpatialPixelsDataFrame for inlabru extract
  bat_spp <- as(bat_proj, "SpatialPixelsDataFrame")
    
  # Extract season and year specific sea ice values at spp coords
  v <- over(spp, bat_spp)
  names(v) <- "depth"
  
  if (any(is.na(v$depth))) {
      v$depth <- inlabru:::bru_fill_missing(bat_spp, spp, v$depth)
    }
  return(v)
}

# ends
