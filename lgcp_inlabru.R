###############################################
### Attempt to fit harp model using inlabru ###
###############################################

# load libraries
require(INLA)
require(inlabru)
require(tidyverse)
require(viridis)
require(patchwork)

# Load a random subsample of ~2500 animal locations
dat <- readRDS("data/harps2500_indexed.rds")

# need to convert this into the right format for inlabru
coordinates(dat) <- cbind(dat$x, dat$y)

prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"
locs <- sp::SpatialPoints(coords = cbind(dat$x, dat$y), proj4string = CRS(prj))

###########################################################################################
### Start simple with a non-barrier lgcp model with no time structure and no covariates ###
###########################################################################################

# simple boundary and mesh
bnd <- inla.nonconvex.hull(points = locs)

# Define the parameters of the boundary mesh
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes ~ less computation time
mesh <- inla.mesh.2d(boundary = bnd,
                     max.edge = c(1, 5) * max.edge,
                     cutoff = 50,
                     offset = c(max.edge, bound.outer),
                     crs = CRS(prj))

# Matern covariance on space
pcmatern <- inla.spde2.pcmatern(mesh,
                                prior.sigma = c(0.1, 0.01),
                                prior.range = c(5, 0.01))

# integration points
ips <- ipoints(mesh)

# model formula - spatial spde and intercept
m1 <- coordinates ~ mySmooth(coordinates, model = pcmatern) + 
                    myDepth(f.bat(x, y), model = "linear") +
                    Intercept(1)

fit1 <- lgcp(m1,
            data = dat,
            samplers = bnd,
            domain = list(coordinates = mesh),
            ips = ips,
            options = list(control.inla = list(int.strategy = "eb")))

summary(fit1)

# Predict the spatial intensity surface
lambda1 <- predict(fit1, pixels(mesh), ~ exp(mySmooth + Intercept))

# Plot the intensity
p1 <- ggplot() +
  gg(lambda1) +
  scale_fill_viridis() +
  coord_fixed()

#############################
### Try the barrier model ###
#############################

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

# create the matern object
barrier.model = inla.barrier.pcmatern(mesh,
                                      barrier.triangles = barrier.triangles,
                                      prior.range = c(50, .1),
                                      prior.sigma = c(10, 0.01))

# integration points
ips <- ipoints(mesh)

# model formula barrier spde and intercept
m2 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh)) + 
                    Intercept(1)

fit2 <- lgcp(m2,
             data = dat,
             samplers = bound,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb")))
summary(fit2)

# Predict the spatial intensity surface
lambda2 <- predict(fit2, pixels(mesh), ~ exp(mySmooth + Intercept))

# Plot the intensity
p2 <- ggplot() +
  gg(lambda2) +
  scale_fill_viridis() +
  coord_fixed()

p1 + p2

deltaIC(fit1, fit2)

####################################################
### Now add the seasonal discrete time component ###
####################################################

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
ips_season <- ipoints(domain = 1:4, name="season")
ips <- cprod(ips_space, ips_season)

# model formula
m3 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model="ar1")) + 
  Intercept(1)

fit3 <- lgcp(m3,
             data = dat,
             samplers = bound,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

summary(fit3)

# compare models
deltaIC(fit1, fit2, fit3)

# generate predictions
ppxl <- pixels(mesh, mask = bound)
ppxl_all <- cprod(ppxl, data.frame(season = seq_len(4)))

lambda3 <- predict(fit3, ppxl_all, ~ data.frame(season = 1:4, lambda = exp(mySmooth + Intercept)))

p3 <- ggplot() +
  gg(lambda3, aes(fill = mean)) +
  scale_fill_viridis() +
  facet_wrap(~season) +
  coord_equal()
p3

##############################
### add the year component ###
##############################
ips_space <- ipoints(sampler = bound, domain = mesh)
ips_season <- ipoints(domain = 1:4, name="season")
ips <- cprod(ips_space, ips_season)
# if this runs ok then it should be easy to add sea ice
ips_year <- ipoints(domain = c(1, 2, 3, 5), name = "year")
ips <- cprod(ips_space, ips_season, ips_year)

# model formula
m4 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model="ar1")) + 
  Intercept(1)

fit4 <- lgcp(m4,
             data = dat,
             samplers = bound,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

summary(fit4)

# need to add year index to prediction object
ppxl_all <- cprod(ppxl, data.frame(season = seq_len(4), year = c(1, 2, 3, 5)))
lambda4 <- predict(fit4, ppxl_all, ~ data.frame(season = 1:4, lambda = exp(mySmooth + Intercept)))

p4 <- ggplot() +
  gg(lambda4, aes(fill = mean)) +
  scale_fill_viridis() +
  facet_wrap(~season) +
  coord_equal()
p4

########################################
### Add linear sea ice concentration ###
########################################

require(raster)
source("R/f.ice_single.R")
ice <- readRDS("data/ice_extend_proj.rds")
ice <- mean(ice)
#ice <- ice - mean(ice, na.rm = T)

values(ice)[values(ice) < 0] <- 0 # force interpolated values less than 0 to be 0

# model formula
# need to check with Finn the difference between Intercept, Intercept(1) and no intercept term?
m5 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1")) +
                    myIce(f.ice(x, y), model = "linear") +
                    Intercept(1)

ice_proj <- ice
ice_mesh <- inla.mesh.1d(loc = seq(from = 0, to = 100, by = 2.5))
ice_spde <- inla.spde2.pcmatern(ice_mesh, prior.range = c(1000, NA), prior.sigma = c(1, NA)) 

#ice_spde <- inla.spde2.matern(ice_mesh, constr = F) 

m6 <- coordinates ~ myIce(f.ice_single(x, y),
                          model = ice_spde) +
                    Intercept(1)

m7 <- coordinates ~ myIce(f.ice_single(x, y),
                          model = "rw2",
                          mapper = bru_mapper(ice_mesh, indexed = F))

m7 <- coordinates ~ myIce(f.ice(x, y),
                          model = "rw2",
                          mapper = bru_mapper(ice_mesh, indexed = F)) +
                    mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model="ar1"))


         
fit8 <- lgcp(m7,
             data = dat,
             samplers = bound,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

lambda6 <- predict(fit8, ppxl_all, data.frame(lambda = exp(mySmooth + myIce + Intercept)))

p1 <- ggplot() + 
  geom_ribbon(aes(x = ice, ymin = q0.025, ymax = q0.975), data = lambda6, fill = "grey") +
  geom_line(aes(x = ice, y = mean), data = lambda6)
p2 <- ggplot() + 
  geom_ribbon(aes(x = ID, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`)), fill = "grey", data = fit7$summary.random$myIce) + 
  geom_line(aes(x = ID, y = exp(mean)), data = fit7$summary.random$myIce) +
  ylim(0, 10)

p1 + p2





summary(fit5)
deltaIC(fit4, fit5)

# need to add year index to prediction object
ppxl_all <- cprod(ppxl, data.frame(season = seq_len(4), year = c(1, 2, 3, 5)))
lambda5 <- predict(fit5, ppxl_all, ~ data.frame(season = 1:4, lambda = exp(mySmooth + myIce + Intercept)))

p5 <- ggplot() +
  gg(lambda6, aes(fill = mean)) +
  scale_fill_viridis() +
  facet_wrap(~season) +
  coord_equal()
p5

########################################
### Non-linear sea ice concentration ###
########################################

#ice_mesh <- inla.mesh.1d(loc = quantile(dat$ice, probs = seq(0, 1, length = 9), na.rm = T))

#ice_mesh <- inla.mesh.1d(loc = c(0, 1, 2.5, 5, 10, 20, 40, 60, 80, 100), boundary = "free")

ice_mesh <- inla.mesh.1d(loc = seq(from = 0, to = 100, by = 1), boundary = "free")
#ice_spde <- inla.spde2.pcmatern(ice_mesh, prior.range = c(1, 0.01), prior.sigma = c(10, 0.01)) 

# model formula
# need to check with Finn the difference between Intercept, Intercept(1) and no intercept term?
m6 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1")) +
  myIce(f.ice(x, y), model = ice_spde) +
  Intercept(1)

m6 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model = "ar1")) +
                    myIce(f.ice(x, y),
                          model = "rw2",
                          mapper = bru_mapper(ice_mesh, indexed = F)) +
                    Intercept(1)

#m6 <- coordinates ~ myIce(f.ice(x, y),
#                          model = "rw2",
#                          mapper = bru_mapper(ice_mesh, indexed = F)) + Intercept(1)

fit6 <- lgcp(m6,
             data = dat,
             samplers = bound,
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))


preds <- predict(fit8, formula = ~ mySmooth + myIce + Intercept)

ggplot() +
  geom_line(aes(x = ice, y = mean), data = pred_ice_spde)

post_range <- spde.posterior(fit6, name = "myIce", what = "range")
plot(post_range)


summary(fit6)
deltaIC(fit4, fit5, fit6)




require(sf)

preds <- rbind(mesh$loc, mesh$loc, mesh$loc, mesh$loc)
preds <- cbind(preds, fit3$summary.random$mySmooth$mean)
preds <- data.frame(preds)

names(preds) <- c("x", "y", "foo", "fit")
preds$season <- rep(1:4, each = mesh$n)
ggplot() + geom_point(aes(x = x, y = y, colour = fit), data = preds) + scale_colour_viridis() + facet_wrap(~season)




m2 <- coordinates ~ mySmooth(coordinates,
                             model = pcmatern,
                             group = season,
                             ngroup = 4,
                             control.group = list(model="ar1")) + season + Intercept(1)

fit2 <- lgcp(m2,
            data = dat,
            samplers = bnd,
            domain = list(coordinates = mesh),
            ips = ips)

summary(fit2)





### Try the barrier model
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

# create the matern object
barrier.model = inla.barrier.pcmatern(mesh,
                                      barrier.triangles = barrier.triangles,
                                      prior.range = c(50, .1),
                                      prior.sigma = c(10, 0.01))

# integration points
ips <- ipoints(mesh)

# model formula
m3 <- coordinates ~ mySmooth(coordinates, model = barrier.model, mapper = bru_mapper(mesh)) + Intercept(1)

fit3 <- lgcp(m3,
            data = dat,
            samplers = bound,
            domain = list(coordinates = mesh),
            ips = ips,
            options = list(control.inla = list(int.strategy = "eb")))
summary(fit3)
            
# Predict the spatial intensity surface
lambda <- predict(fit3, pixels(mesh), ~ exp(mySmooth + Intercept))

# Plot the intensity
ggplot() +
  gg(lambda) +
  scale_fill_viridis() +
  coord_fixed()






### try ice with the f_ice function..
source("R/f.ice.R")
ice <- readRDS("data/ice_proj.rds")
ice <- ice - mean(ice, na.rm = TRUE)



# simple boundary and mesh
bnd <- inla.nonconvex.hull(points = locs)

# Define the parameters of the boundary mesh
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes ~ less computation time
mesh <- inla.mesh.2d(boundary = bnd,
                     max.edge = c(1, 5) * max.edge,
                     cutoff = 50,
                     offset = c(max.edge, bound.outer),
                     crs = CRS(prj))

# Matern covariance on space
pcmatern <- inla.spde2.pcmatern(mesh,
                                prior.sigma = c(0.1, 0.01),
                                prior.range = c(5, 0.01))

# integration points
ips <- ipoints(mesh)

# model formula
m4 <- coordinates ~ myIce(f.ice(x, y), model = ice_spde) + 
                    mySmooth(coordinates, model = pcmatern) +
                    Intercept

fit4 <- lgcp(m4,
             data = dat,
             samplers = bound,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))
summary(fit4)

lambda4 <- predict(fit4, pixels(mesh), ~ exp(mySmooth + myIce + Intercept))

lambda5 <- predict(fit4, pixels(ice_mesh), ~ exp(myIce + Intercept))

require(patchwork)
require(viridis)
p1 <- ggplot() +
  gg(lambda) +
  scale_fill_viridis() +
  coord_fixed()
p2 <- ggplot() +
  gg(lambda4) +
  scale_fill_viridis() +
  coord_fixed()

p1 + p2


deltaIC(fit, fit4)


plot(seq(from = -65, to = 65, by = 10), exp(fit4$summary.random$myIce$mean), ylim = c(0, 2))
lines(seq(from = -65, to = 65, by = 10), exp(fit4$summary.random$myIce$`0.025quant`))
lines(seq(from = -65, to = 65, by = 10), exp(fit4$summary.random$myIce$`0.975quant`))




emdl <- coordinates ~ mySmooth(coordinates, model = matern) + Intercept

fit <- lgcp(emdl, dat_sp, samplers = bound, domain = list(coordinates = mesh))




require(inlabru)


ggplot() +
  gg(ice) +
  gg(boundary, alpha = 0) +
  coord_fixed()



ggplot() + 
  geom_sf(aes(fill = w), data = dmesh %>% st_as_sf()) +
  geom_sf(aes(), fill = NA, colour = "white", data = bound %>% st_as_sf()) +
  scale_fill_viridis() +
  theme_bw()


###

bathy <- readRDS("data/bathy.rds")
# set land to zero


prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"
bat_proj <- raster::projectRaster(bathy, crs = CRS(prj), res = 25)
bat_proj[bat_proj >= 0] <- NA
source("R/f.bat.R")

bat_proj <- bat_proj * -1
bat_proj <- log(bat_proj)

require(INLA)
require(inlabru)
require(tidyverse)
require(viridis)
require(patchwork)

# Load a random subsample of ~2500 animal locations
dat <- readRDS("data/harps2500_indexed.rds")

# need to convert this into the right format for inlabru
coordinates(dat) <- cbind(dat$x, dat$y)

prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"
locs <- sp::SpatialPoints(coords = cbind(dat$x, dat$y), proj4string = CRS(prj))

###########################################################################################
### Start simple with a non-barrier lgcp model with no time structure and no covariates ###
###########################################################################################

# simple boundary and mesh
bnd <- inla.nonconvex.hull(points = locs)

# Define the parameters of the boundary mesh
max.edge = max(c(diff(range(dat$x)), diff(range(dat$y))))/15
bound.outer = 1000

# Changing cutoff is an effective way to change number of mesh nodes
# Fewer mesh nodes ~ less computation time
mesh <- inla.mesh.2d(boundary = bnd,
                     max.edge = c(1, 5) * max.edge,
                     cutoff = 50,
                     offset = c(max.edge, bound.outer),
                     crs = CRS(prj))

# Matern covariance on space
pcmatern <- inla.spde2.pcmatern(mesh,
                                prior.sigma = c(0.1, 0.01),
                                prior.range = c(5, 0.01))

# integration points
ips <- ipoints(mesh)

# model formula - spatial spde and intercept
m1 <- coordinates ~ mySmooth(coordinates, model = pcmatern) + 
                    myDepth(f.bat(x, y), model = "linear") +
                    Intercept(1)

fit1 <- lgcp(m1,
             data = dat,
             samplers = bnd,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

summary(fit1)

# Predict the spatial intensity surface
lambda1 <- predict(fit1, pixels(mesh), ~ exp(mySmooth + Intercept))

# Plot the intensity
p1 <- ggplot() +
  gg(lambda1) +
  scale_fill_viridis() +
  coord_fixed()

bat_mesh <- inla.mesh.1d(loc = seq(from = -6, to = 8, by = 1), boundary = "free")
bat_spde <- inla.spde2.pcmatern(bat_mesh, prior.range = c(5, 0.01), prior.sigma = c(10, 0.01)) 

#using log10 transform and range of 1 and sigma of 3 got a wiggly line...
# can go back to that if this doesn't work
# easier to back transform from log using exp than if transformation is log10...
# difficulty is scale of projected bathymetry
# the raw raster doesn't go down as deep as the projected raster - must be an issue with the interpolation
# check to make sure the model fit isn't estimating the marianna trench!!!


# model formula barrier spde and intercept
m2 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh)) + 
                    myDepth(f.bat(x, y),
                            model = "rw2",
                            mapper = bru_mapper(bat_mesh,
                                                indexed = F)) +
                    Intercept(1)

m2 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh)) + 
                    myDepth(f.bat(x, y),
                            model = bat_spde) +
  Intercept(1)


fit2 <- lgcp(m2,
             data = dat,
             samplers = bnd,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

plot(exp(bat_mesh$mid), exp(fit2$summary.fixed$mean + (fit2$summary.random$myDepth$mean*bat_mesh$mid)), type = "l", xlim = c(0, 10000))

lambda1 <- predict(fit2, pixels(mesh), ~ exp(mySmooth + myDepth + Intercept))
lambda2 <- predict(fit2, pixels(mesh), ~ exp(mySmooth + Intercept))

lambda1 <- predict(fit2, pixels(mesh), ~ exp(mySmooth + Intercept))


bat_proj <- bat_proj2


### try both bathy and static ice

bat_mesh <- inla.mesh.1d(loc = seq(from = -6, to = 8, by = 1), boundary = "free")
bat_spde <- inla.spde2.pcmatern(bat_mesh, prior.range = c(5, 0.01), prior.sigma = c(10, 0.01)) 

source("R/f.ice_single.R")
ice <- readRDS("data/ice_proj.rds")
ice <- raster::mean(ice)
ice_proj <- ice
ice_proj[ice_proj < 0] <- 0

ice_mesh <- inla.mesh.1d(loc = seq(from = 0, to = 100, by = 1), boundary = "free")
ice_spde <- inla.spde2.pcmatern(ice_mesh, prior.range = c(25, 0.01), prior.sigma = c(50, 0.01)) 

m3 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh)) + 
                    myDepth(f.bat(x, y),
                            model = bat_spde) +
                    myIce(f.ice_single(x, y),
                          model = ice_spde) +
                    Intercept(1)

fit4 <- lgcp(m3,
             data = dat,
             samplers = bnd,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

ice <- readRDS("data/ice_proj.rds")
ice_proj <- ice
ice_proj[ice_proj < 0] <- 0

m4 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh)) + 
                    myDepth(f.bat(x, y),
                            model = bat_spde) +
                    myIce(f.ice(x, y),
                          model = ice_spde) +
                    Intercept(1)

fit5 <- lgcp(m4,
             data = dat,
             samplers = bnd,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))



#################################################
### non linear bathy, seasonal ice and AR1... ###
#################################################
# this will take 10 hours and the fit is terrible

# could try starting values:
bru_initial = c(1.8722477, 8.2904219, 0.8382415, 58.1938352, 2.6262642, 174.5949562, 3.3875526)


ice <- readRDS("data/ice_proj.rds")
ice <- raster::stackApply(ice, indices = seq(1:4), fun = mean, na.rm = T)
ice_proj <- ice
ice_proj[ice_proj < 0] <- 0
#f_ice helper function limited to seasonal average...
source("R/f.ice_seasonal.R")

# bathymetry
bathy <- readRDS("data/bathy.rds")
prj = "+proj=laea +lat_0=75 +lon_0=-25 +x_0=0 +y_0=0 +units=km +no_defs +ellps=WGS84"
bat_proj <- raster::projectRaster(bathy, crs = CRS(prj), res = 25)
bat_proj[bat_proj >= 0] <- NA
source("R/f.bat.R")
bat_proj <- bat_proj * -1
bat_proj <- log(bat_proj)

# integration points
ips_space <- ipoints(sampler = bound, domain = mesh)
ips_season <- ipoints(domain = 1:4, name="season")
ips <- cprod(ips_space, ips_season)
# if this runs ok then it should be easy to add sea ice
#ips_year <- ipoints(domain = c(1, 2, 3, 5), name = "year")
#ips <- cprod(ips_space, ips_season, ips_year)

ice_mesh <- inla.mesh.1d(loc = seq(from = 0, to = 100, by = 1), boundary = "free")
ice_spde <- inla.spde2.pcmatern(ice_mesh, prior.range = c(25, 0.01), prior.sigma = c(25, 0.01)) 

bat_mesh <- inla.mesh.1d(loc = seq(from = -6, to = 8, by = 1), boundary = "free")
bat_spde <- inla.spde2.pcmatern(bat_mesh, prior.range = c(5, 0.01), prior.sigma = c(10, 0.01)) 

m5 <- coordinates ~ mySmooth(coordinates,
                             model = barrier.model,
                             mapper = bru_mapper(mesh),
                             group = season,
                             ngroup = 4,
                             control.group = list(model="ar1")) + 
                    myDepth(f.bat(x, y),
                            model = bat_spde) +
                    myIce(f.ice(x, y),
                          model = ice_spde) +
                    Intercept(1)

fit6 <- lgcp(m5,
             data = dat,
             samplers = bnd,
             domain = list(coordinates = mesh),
             ips = ips,
             options = list(control.inla = list(int.strategy = "eb"),
                            verbose = T))

ppxl <- pixels(mesh, mask = bound)
ppxl_all <- cprod(ppxl, data.frame(season = seq_len(4)))
lambda6 <- predict(fit6, ppxl_all, ~ data.frame(season = 1:4, lambda = exp(mySmooth + myDepth + myIce + Intercept)))
###

# this complicated model estimated the same result with ice... so should be able to tweak simpler model ice...



