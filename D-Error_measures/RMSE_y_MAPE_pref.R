#*********************************
# Error measures                 #
# RMSE y MAPE                    #
# CPUE indices                   #
# Preferential sampling          #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Packages ---------------------------
library(Metrics)
library(MLmetrics)

# Median for the simulated biomass and effort each year ---------------------------
load("./A-Simulation/simulated_data.RData")
biomasa_simulada <- list()
for (i in 1:k) {
  biomasa_simulada[[i]] <- median(data_variables[[i]]$biomasa)
} 

biomasa_simulada <- unlist(biomasa_simulada)

esfuerzo_por_anyo <- list()
for (i in 1:k) {
  esfuerzo_por_anyo[[i]] <- median(data_pref_final$esfuerzo[data_pref_final$tiempo == i])
}
esfuerzo_por_anyo <- unlist(esfuerzo_por_anyo)

# Clear the index to get the biomass (paper formulation) ---------------------------
load("./C-Estimated_series/CPUE_series.RData")
glm_f_bio <- (glm_f_CPUE * esfuerzo_por_anyo)/q_pref # frequentist GLM
gam_f_bio <- (gam_f_CPUE * esfuerzo_por_anyo)/q_pref # frequentist GAM 
glm_b_bio <- (glm_b_CPUE * esfuerzo_por_anyo)/q_pref # Bayesian GLM 
gam_b_bio <- (gam_b_CPUE * esfuerzo_por_anyo)/q_pref # Bayesian GAM 
geo_b_bio <- (geo_b_CPUE * esfuerzo_por_anyo)/q_pref # Geostatistical model
geo_pref_b_bio <- (geo_pref_b_CPUE * esfuerzo_por_anyo)/q_pref # Marked point pattern

# RMSE ---------------------------
rmse(glm_f_bio, biomasa_simulada)
rmse(gam_f_bio, biomasa_simulada)
rmse(glm_b_bio, biomasa_simulada)
rmse(gam_b_bio, biomasa_simulada)
rmse(geo_b_bio, biomasa_simulada)
rmse(geo_pref_b_bio, biomasa_simulada)

# MAPE ---------------------------
MAPE(glm_f_bio, biomasa_simulada)
MAPE(gam_f_bio, biomasa_simulada)
MAPE(glm_b_bio, biomasa_simulada)
MAPE(gam_b_bio, biomasa_simulada)
MAPE(geo_b_bio, biomasa_simulada)
MAPE(geo_pref_b_bio, biomasa_simulada)
