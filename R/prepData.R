##########################################
### Prepare the data for the INLA LGCP ###
##########################################

######################
### Load libraries ###
######################
require(INLA)
require(inlabru)
require(tidyverse)
require(sf)
require(viridis)
require(raster)
require(mapr)
require(rgeos)
source("R/spde-book-functions.R")

################################
### Format the location data ###
################################

# Here is a random subsample of 2500 animal locations
dat <- readRDS("data/harps2500_indexed.rds")

# Append identity of which harp seal population animals come from
# Load in the population reference table
pop_ref <- readRDS("data/populations.rds")

# Match the ids with pop ref and take the corresponding population reference
dat <- dat %>% left_join(pop_ref, by = c("id" = "ref"))
dat <- dat %>% dplyr::select("id", "date", "lon", "lat", "x", "y", "year", "month", "index", "year_i", "season", "location")

# forgot to append all seals to the population table...
dat <- dat %>% mutate(population = case_when(location == "Newfoundland" ~ 1,
                                             location == "West Ice" ~ 2,
                                             location == "East Ice" ~ 3,
                                             id == "hp5-L764-18" ~ 3,
                                             id == "hp5-L766-18" ~ 3,
                                             TRUE ~ 1)) # hp6 ~ 1

#################################################
### Define the spatial mesh and barrier model ###
#################################################

# Define an Albers equal area projection using kilometres as units
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"

# Generate a land shapefile for generating a map later on
land <- mapr::mapr(dat,
                   prj,
                   buff = 2000)

# Generate simple boundary for the inla mesh
b <- mapr::meshr(dat,
                 prj,
                 buff = 500,
                 keep = 0.5,
                 Neumann = F)

# Remove any animal locations that fall on land
in.water = over(b, SpatialPoints(cbind(dat$x, dat$y), proj4string = CRS(prj)), returnList=T)[[1]]
dat <- dat[in.water,]

# Define the parameters of the boundary mesh
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes ~ less computation time
mesh <- inla.mesh.2d(boundary = b,
                     max.edge = c(1, 5) * max.edge,
                     cutoff = 50,
                     offset = c(max.edge, bound.outer),
                     crs = CRS(prj))

# Construct Voroni polygons
dmesh <- book.mesh.dual(mesh)
proj4string(dmesh) <- prj

require(rgeos)
b <- gBuffer(b, byid=TRUE, width=0) # use gBuffer to stop b is invalid errors

# Save objects
saveRDS(land, "data/land_sf.rds")
saveRDS(b, "data/boundary.rds")
saveRDS(mesh, "data/mesh.rds")
saveRDS(dmesh, "data/dmesh.rds")

###########################################
### Append space-time sea ice covariate ###
###########################################

# Load and stack all 20 raster layers together
ice <- raster::stack(raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec94_Nov99"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec99_Nov04"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec04_Nov09"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec09_Nov14"),
                     raster::stack("NSIDC Sea Ice/Seasonal_NSIDC_Dec14_Nov19"))

# Extract covariate to mesh...
# should this be averaged across the voroni polygon or just extracted to mesh node?
# could use dmesh instead...?
# need to do this for each year and each season...
mesh_ice <- extract(ice,
                    dmesh %>%
                      st_as_sf() %>%
                      st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                    fun = mean, na.rm = T)

# Repeat for data points
# Extract space-time matched data using unique times index
dat$ice <- NA
for (i in sort(unique(dat$index))){
  dat$ice[dat$index == i] <- raster::extract(raster::subset(ice, i),
                                             dat[dat$index == i,] %>% st_as_sf(coords = c("x", "y"), crs = prj) %>%
                                               st_transform("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"),
                                             method = "bilinear")
}

saveRDS(dat, "data/harps2500_indexed.rds")
saveRDS(mesh_ice, "data/mesh_ice.rds")

# ends