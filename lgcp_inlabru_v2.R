##################################################################
### Harp seal inlabru lgcp with heirarchical sea ice covariate ###
##################################################################

# 2022-04-25

# An example script to run an inlabru lgcp model on the harp seal telemetry data
# Model setup is space-time lgcp
# Barrier model
# Seasonal effect on distribution
# Multiple years of data added to ips seperately due to differences in sea ice characteristics
# Heirarchically centred sea ice covariate - seperating out spatial pattern, seasonal pattern and decline in ice over time

#####################################
### Load libraries and point data ###
#####################################

# libraries
require(INLA)
require(inlabru)
require(tidyverse)
require(viridis)
require(patchwork)
require(sf)
require(raster)

# Load a random subsample of ~2500 animal locations
dat <- readRDS("data/harps2500_indexed.rds")

# need to convert this into the right format for inlabru
coordinates(dat) <- ~x+y

# define projection
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"

############################
### Set up barrier model ###
############################

bound <- readRDS("data/boundary.rds")
mesh <- readRDS("data/mesh.rds")

# Set up boundary matern model
tl = length(mesh$graph$tv[,1]) # the number of triangles in the mesh

# Compute triangle positions
posTri = matrix(0, tl, 2) 
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri, proj4string = CRS(prj))

normal = unlist(over(bound, posTri, returnList = T)) # check which mesh triangles are inside the normal area
barrier.triangles = setdiff(1:tl, normal)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)

# create the matérn object
barrier.model = inla.barrier.pcmatern(mesh,
                                      barrier.triangles = barrier.triangles,
                                      prior.range = c(50, .1),
                                      prior.sigma = c(10, 0.01))

###########################################################################
### Now add the seasonal discrete time component and the year component ###
###########################################################################

# this highlighted an issue with ips integration
# if mesh is projected and covers the pole this will cause an issue with internal reprojection

# one solution is to - 
# replace ips with the spatial points of the mesh
# add weights according to the dmesh integration weights from the voroni polygons
#w <- readRDS("data/voroni_weights.rds")
#ips_space <- SpatialPointsDataFrame(mesh$loc, data = data.frame(weight = w))

# alternatively strip the projection from the boundary and mesh files
# this will force INLA to work with them in the native (correct) projection
proj4string(bound) <- sp::CRS(NA_character_)
mesh$crs <- sp::CRS(NA_character_)

# then use these to create a space-time indexed ips object
ips_space <- ipoints(sampler = bound, domain = mesh)
ips_season <- ipoints(domain = 1:4, name = "season")
ips_year <- ipoints(domain = c(1, 2, 3, 5), name = "year_i")
ips <- cprod(ips_space, ips_season, ips_year)

# adjust the exposure of the ips object based on tagging effort
source("R/adjust_exposure.R")
ips <- adjust_exposure(dat, ips)

#############################################################
### Put this all together in an intercept only lgcp model ###
#############################################################

# set up a prior on the AR1 process
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))

# model formula
m0 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1",
                                                  hyper = list(theta = pcrho))) +
  Intercept(1)

fit0 <- lgcp(m0,
             data = dat,
             samplers = bound,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            control.mode = list(theta = c(.7, 7, .7), restart = TRUE),
                            verbose = T))
summary(fit0)

# need to add year index to prediction object
ppxl <- pixels(mesh, mask = bound)
ppxl_all <- cprod(ppxl, data.frame(season = seq_len(4), year_i = c(1, 2, 3, 5)))
lambda4 <- predict(fit4, ppxl_all, ~ data.frame(season = season, lambda = exp(mySmooth + Intercept)))

p1 <- ggplot() +
  gg(lambda4, aes(fill = mean)) +
  scale_fill_viridis() +
  facet_wrap(~season) +
  coord_equal()
p1

##########################################
### Add hierarchical sea ice covariate ###
##########################################

# sea ice evaluator functions
source("R/f.ice_single.R")
source("R/f.ice_seasonal.R")
source("R/f.ice_residual.R")

# load hierarchical ice components
x_bar <- readRDS("data/x_bar.rds") # average ice over all time - scalar
x_bar_sc <- readRDS("data/x_bar_sc.rds") # average across all time - single spatial field
x_bar_smc <- readRDS("data/x_bar_smc.rds") # average across each season - 4 spatial fields
x_smy_resid <- readRDS("data/x_smy_resid.rds") # residual space-time deviation - year effect

# project to same projection as data
ice_proj <- projectRaster(x_bar_sc, crs = prj)
ice <- projectRaster(x_bar_smc, crs = prj)
ice_resid <- projectRaster(x_smy_resid, crs = prj)

# set up model formulas
# test whether distribution is influenced by long-term trend in ice (the space-time residual)
m1 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1",
                                                  hyper = list(theta = pcrho))) +
  myIceSingle(f.ice_single(x, y), model = "linear") + 
  myIceSeasonal(f.ice_seasonal(x, y, season), model = "linear") + 
  Intercept(1)

m2 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1",
                                                  hyper = list(theta = pcrho))) +
  myIceSingle(f.ice_single(x, y), model = "linear") + 
  myIceSeasonal(f.ice_seasonal(x, y, season), model = "linear") + 
  myIceResidual(f.ice_residual(x, y, season, year_i), model = "linear") + 
  Intercept(1)

# fit models
fit1 <- lgcp(m1,
             data = dat,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            control.mode = list(theta = c(.7, 7, .7), restart = TRUE),
                            verbose = T))
fit2 <- lgcp(m2,
             data = dat,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            control.mode = list(theta = c(.7, 7, .7), restart = TRUE),
                            verbose = T))

# check outputs
summary(fit1)
summary(fit2)

deltaIC(fit1, fit2)

# ends