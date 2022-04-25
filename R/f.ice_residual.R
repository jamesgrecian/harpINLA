


# function to extract sea ice raster data to xy coordinates at time point t
# sea ice input is a raster stack with 20 layers - 4 season * 5 years
# function input is season (1:4) and year (1:5) from which time point t is calculated (1:20)

f.ice_residual <- function(x, y, season, year) {
  
  # convert season and year to indices for raster stack
  t <- (year - 1) * 4 + season
  
  # extract space time matched sea ice to points
  # loop through unique layers
  v <- rep(NA, length(t))
  
  for(tt in unique(t)){
    
    # turn coordinates into SpatialPoints object:
    # with the appropriate coordinate reference system (CRS)
    # and the appropriate time indices
    spp <- SpatialPoints(coords = cbind(x = x[t == tt], y = y[t == tt]), proj4string = fm_sp_get_crs(ice_resid))
    
    # Subset raster stack at correct time index
    ice_resid_sub <- raster::subset(ice_resid, tt)
    # Convert to SpatialPixelsDataFrame for inlabru extract
    ice_resid_sub <- as(ice_resid_sub, "SpatialPixelsDataFrame")
    
    # Extract season and year specific sea ice values at spp coords
    result <- over(spp, ice_resid_sub)
    names(result) <- "ice"
    
    if (any(is.na(result$ice))) {
      result$ice <- inlabru:::bru_fill_missing(ice_resid_sub, spp, result$ice)
    }
    
    v[t == tt] <- result$ice
  }
  return(v)
}

# ends

