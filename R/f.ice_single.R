# function to extract bathymetry data to xy coordinates

f.ice_single <- function(x, y) {
  
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(coords = cbind(x = x, y = y), proj4string = fm_sp_get_crs(ice_proj))
  
  # Convert to SpatialPixelsDataFrame for inlabru extract
  ice_spp <- as(ice_proj, "SpatialPixelsDataFrame")
  
  # Extract season and year specific sea ice values at spp coords
  v <- over(spp, ice_spp)
  names(v) <- "ice"
  
  if (any(is.na(v$ice))) {
    v$ice <- inlabru:::bru_fill_missing(ice_spp, spp, v$ice)
  }
  return(v)
}

# ends
