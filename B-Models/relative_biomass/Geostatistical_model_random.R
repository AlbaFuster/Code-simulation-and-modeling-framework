#*********************************
# Geostatistical model           #
# Relative biomass data          #
# Random sampling                #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Load data ---------------------------
load("./A-Simulation/simulated_data.RData")

# Packages ---------------------------
library(INLA)
library(raster)

# Eliminate the values of the simulated biomass ---------------------------
data_random_final <- data_random_final[, -1] 

# Mesh ---------------------------
loc <- cbind(data_random_final$xcoord,
             data_random_final$ycoord
             ) # Locations

mesh <- inla.mesh.2d(loc.domain = cbind(c(0,10,10,0,0),c(0,0,10,10,0)), # Limit
                     max.edge = c(0.63, 2.5), # Mesh parameters
                     cutoff = 0.05
)

# Data for the prediction ---------------------------
data_bat <- data.frame(x = loc_xy[, 1],
                       y = loc_xy[, 2],
                       z = variables$x_bat1
                       )
ras_bat <- rasterFromXYZ(data_bat) # Transform bathymetry into a raster
bat_pred <- raster::extract(ras_bat,
                            mesh$loc[, c(1, 2)]
                            ) # Extract values of the mesh

data_pred <- data.frame(batimetria = bat_pred,
                        xcoord = mesh$loc[, 1],
                        ycoord = mesh$loc[, 2]) # Data

data_pred_final <- rbind(
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred,
  data_pred
) # repeat the data k (years) times

data_pred_final$tiempo <- rep(1:15, each = length(data_pred$batimetria)) # add years
data_pred_final$bio_relativa <- rep(NA, length(data_pred_final$batimetria)) # add relative biomass

# INLA.GROUP ---------------------------
data_pred_final$batimetria[which(data_pred_final$batimetria > max(data_random_final$batimetria))] <- NA # Same limits for pred and est
data_pred_final$batimetria[which(data_pred_final$batimetria < min(data_random_final$batimetria))] <- NA

all <- rbind(data_random_final, data_pred_final) # Joint pred and est
igroupbath <- inla.group(all[, 1], n = 20, idx.only = TRUE) # Funcion inla.group
groupbath <- inla.group(all[, 1], n = 20) # Funcion inla.group

allin <- cbind(all, igroupbath, groupbath) # Joint

data_random_final <- cbind(data_random_final,
                           allin$groupbath[1:length(data_random_final$batimetria)],
                           allin$igroupbath[1:length(data_random_final$batimetria)]
                           ) # Add columns inla.group

data_pred_final <- cbind(data_pred_final,
                         allin$groupbath[(length(data_random_final$batimetria) + 1):length(allin$batimetria)],
                         allin$igroupbath[(length(data_random_final$batimetria) + 1):length(allin$batimetria)]
                         ) # Add columns inla.group

names(data_random_final) <- c("batimetria",
                              "xcoord",
                              "ycoord",
                              "bio_relativa",
                              "tiempo",
                              "IgBath",
                              "Igroup"
                              )
names(data_pred_final) <- c("batimetria",
                            "xcoord", 
                            "ycoord", 
                            "tiempo", 
                            "bio_relativa", 
                            "IgBath", 
                            "Igroup"
                            )

# SPDE ---------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(3, 0.5), # P(range < 3) = 0.5
                            prior.sigma = c(1, 0.01) # P(sigma > 1) = 0.01
                            ) 

# Matrix of weights (spatial field) ---------------------------
iset <- inla.spde.make.index("i",
                             n.spde = spde$n.spde,
                             n.group = k
                             )

# Projection matrix ---------------------------
A_est <- inla.spde.make.A(mesh = mesh,
                          loc = cbind(data_random_final$xcoord, data_random_final$ycoord),
                                group = data_random_final$tiempo
                          ) # estimation

A_pred <- inla.spde.make.A(mesh = mesh,
                          loc = cbind(data_pred_final$xcoord, data_pred_final$ycoord),
                          group = data_pred_final$tiempo
) # prediction

# Stack ---------------------------
## Estimation ---------------------------
stack_est <- inla.stack(data = list(y = data_random_final$bio_relativa, link = 1), 
                        A = list(A_est, 1, 1, 1), 
                        effects = list(iset,
                          batimetria = data_random_final$IgBath,
                          tiempo = data_random_final$tiempo,
                          intercept = rep(1, length(data_random_final$bio_relativa))
                        ),
                        tag = "est"
                        )

## Prediction ---------------------------
stack_pred <- inla.stack(data = list(y = data_pred_final$bio_relativa, link = 1), # response variable NA
                         A = list(A_pred, 1, 1, 1), 
                         effects = list(iset,
                           batimetria = data_pred_final$IgBath,
                           tiempo = data_pred_final$tiempo,
                           intercept = rep(1, length(data_pred_final$bio_relativa))
                         ),
                         tag = "pred"
)

## Joint both stacks ---------------------------
stack <- inla.stack(stack_est, stack_pred)

# PC-prior rho ---------------------------
h_spec <- list(rho = list(prior = "pc.cor1", param = c(0, 0.9))) # P(cor > 0 = 0.9)

# Formula ---------------------------
formula <- y ~ -1 +
               intercept +
               f(batimetria, model = "rw2") + # smooth function for the bathyemtry
               f(i, model = spde, group = i.group, # spatio-temporal effect
                 control.group = list(model = "ar1", hyper = h_spec)) +
               f(tiempo, model = "rw1") # temporal effect

# Fit the model ---------------------------
set.seed(12346)
modelo <- inla(formula,
               data = inla.stack.data(stack),
               family = "gamma",
               control.predictor = list(
                 compute = TRUE,
                 A = inla.stack.A(stack), link = 1),
               verbose = TRUE,
               control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE),
               num.threads = 2
               )

# Prediction ---------------------------
index.p <- inla.stack.index(stack, tag = "pred")$data

## Median of the posterior predictive distribution ---------------------------
output_median <- exp(modelo$summary.linear.predictor[index.p,"0.5quant"])
bio_median <- output_median

prj <- inla.mesh.projector(mesh,
  xlim = range(loc_xy[, 1]),
  ylim = range(loc_xy[, 2])
)

k <- 15
m_mesh <- mesh$n

m_prj_median <- lapply(1:k, function(j) {
  r_median <- inla.mesh.project(
    prj,
    bio_median[1:m_mesh + (j - 1) * m_mesh]
  )
  return(r_median)
})

## Sd of the posterior predictive distribution ---------------------------
output_sd <- exp(modelo$summary.linear.predictor[index.p, "sd"])
bio_sd <- output_sd

m_prj_sd <- lapply(1:k, function(j) {
  r_sd <- inla.mesh.project(
    prj,
    bio_sd[1:m_mesh + (j - 1) * m_mesh]
  )
  return(r_sd)
})

## q0.025 of the posterior predictive distribution ---------------------------
output_25 <- exp(modelo$summary.linear.predictor[index.p, "0.025quant"])
bio_25 <- output_25

m_prj_25 <- lapply(1:k, function(j) {
  r_25 <- inla.mesh.project(
    prj,
    bio_25[1:m_mesh + (j - 1) * m_mesh]
  )
  return(r_25)
})

## q0.975 of the posterior predictive distribution ---------------------------
output_975 <- exp(modelo$summary.linear.predictor[index.p, "0.975quant"])
bio_975 <- output_975

m_prj_975 <- lapply(1:k, function(j) {
  r_975 <- inla.mesh.project(
    prj,
    bio_975[1:m_mesh + (j - 1) * m_mesh]
  )
  return(r_975)
})

# Save models ---------------------------
save.image("./B-Modelos/relative_biomass/random_geostatistic.RData")
