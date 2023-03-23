#*********************************
# Marked point process	         #
# CPUE data                      #
# Preferential sampling          #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Load data ---------------------------
load("./A-Simulation/simulated_data.RData")

# Packages ---------------------------
library(INLA)
library(ggplot2)
library(inlabru)
library(devtools)
library(INLAutils)
library(rgdal)
library(tidyverse)
bru_options_set(inla.mode = "experimental")

# Modify data ---------------------------
data <- data_pref_final
data$batimetria <- as.numeric(scale(data$batimetria))
spdf <- SpatialPointsDataFrame(cbind(data$xcoord,data$ycoord),
                               data[,-c(1,3,4,5,6)]
                               )

# Limit ---------------------------
poly1 <- Polygon(cbind(c(0,10,10,0,0),c(0,0,10,10,0)))
primer_poly <- Polygons(list(poly1), ID = "A")
boundary <- SpatialPolygons(list(primer_poly))

# Mesh ---------------------------
mesh <- inla.mesh.2d(loc.domain =  cbind(c(0,10,10,0,0),c(0,0,10,10,0)),
                                   max.edge = c(0.63, 2.5), 
                                   cutoff = 0.05
                                   )

# Integration points ---------------------------
ips <- ipoints(mesh, group = "tiempo")

# SPDE ---------------------------
spde <- inla.spde2.pcmatern(mesh,
                            prior.sigma = c(1, 0.01), # P(range < 3) = 0.5
                            prior.range = c(4, 0.9) # P(sigma > 1) = 0.01
                            )

# PC-prior rho ---------------------------
prec.prior <- list(prec = list(param = c(0.001, 0.001)))
h_spec <- list(rho = list(prior = "pc.cor1", param = c(0, 0.9))) # P(cor > 0 = 0.9)

# Components ---------------------------
cmp = ~ 
  # Correlated spatial effecto
  spatial(coordinates,model = spde,
                group = tiempo,
                ngroup = k,
                control.group = list(model = "ar1", hyper = h_spec)) +
  # Copy the correlated spatial effect in the Point pattern with a sclae parameter (fixed = FALSE)
  spatialCopy(coordinates, 
              copy = "spatial",
              group = tiempo,
              ngroup = k,
              control.group = list(model = "ar1", hyper = h_spec),
              fixed = FALSE) +
  # Intercept for the point pattern
  lgcpIntercept(1) +
  # Intercept for the CPUE
  Intercept(1) +
  # Temporal trend 
  temporal(tiempo, model = "rw1", hyper = prec.prior) 

lik1 <- like(data = spdf, # likelihood for the CPUE
             family = "gamma",
             formula = CPUE ~ spatial + Intercept + temporal)
lik2 <- like(data = spdf, # likelihood for the point pattern
             family = "cp",
             ips = ips,
             domain = list(coordinates = mesh),
             formula = coordinates + tiempo ~ spatialCopy + lgcpIntercept)

# Fit model ---------------------------
fit <- bru(components = cmp,
           lik1,
           lik2,
           options = list(verbose = TRUE))

# Prediction ---------------------------
ppxl <- pixels(mesh, mask = boundary)
ppxl_all <- cprod(ppxl, data.frame(tiempo = seq_len(k)))
mark <- predict(fit,
                ppxl_all,
                ~ data.frame(tiempo = tiempo,
                             mark = exp(spatial + Intercept + temporal)
                             )
                )

# Save model ---------------------------
save.image("./B-Modelos/marked_point_process.RData")
