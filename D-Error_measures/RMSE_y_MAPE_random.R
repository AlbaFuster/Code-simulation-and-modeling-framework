#*********************************
# Error measures                 #
# RMSE y MAPE                    #
# Relative biomass indices       #
# Random sampling                #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Packages ---------------------------
library(Metrics)
library(MLmetrics)

# Simulated biomass ---------------------------
load("./A-Simulation/simulated_data.RData")
biomasa_simulada <- list()
for (i in 1:k) {
  biomasa_simulada[[i]] <- median(data_variables[[i]]$biomasa)
} 

biomasa_simulada <- unlist(biomasa_simulada)

# Clear the index to get the biomass (paper formulation) ---------------------------
load("./C-Estimated_series/r_biomass_series.RData")
glm_f_bio <- glm_f_bio_r/q_random # frequentist GLM 
gam_f_bio <- gam_f_bio_r/q_random # frequentist GAM 
glm_b_bio <- glm_b_bio_r/q_random # Bayesian GLM 
gam_b_bio <- gam_b_bio_r/q_random # Bayesian GAM 
geo_b_bio <- geo_b_bio_r/q_random # Geostatistical model

# RMSE ---------------------------
rmse(glm_f_bio, biomasa_simulada)
rmse(gam_f_bio, biomasa_simulada)
rmse(glm_b_bio, biomasa_simulada)
rmse(gam_b_bio, biomasa_simulada)
rmse(geo_b_bio, biomasa_simulada)

# MAPE ---------------------------
MAPE(glm_f_bio, biomasa_simulada)
MAPE(gam_f_bio, biomasa_simulada)
MAPE(glm_b_bio, biomasa_simulada)
MAPE(gam_b_bio, biomasa_simulada)
MAPE(geo_b_bio, biomasa_simulada)