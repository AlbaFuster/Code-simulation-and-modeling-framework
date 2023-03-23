#*********************************
# Estimated series               #
# Relative biomass indices       #
# CPUE indices                   #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Relative biomass indices ---------------------------
load("./B-Models/relative_biomass/GLM_GAM_random.RData")
## Frequentist ---------------------------
### GLM ---------------------------
glm_f_bio_r <- list()
for (i in 1:k) {
  glm_f_bio_r[[i]] <- median(data_pred_glm$respuesta[data_pred_glm$tiempo == i])
}
glm_f_bio_r <- unlist(glm_f_bio_r)

### GAM ---------------------------
gam_f_bio_r <- list()
for (i in 1:k) {
  gam_f_bio_r[[i]] <- median(data_pred_gam$respuesta[data_pred_gam$tiempo == i])
}
gam_f_bio_r <- unlist(gam_f_bio_r)

## Bayesian ---------------------------
### GLM --------------------------- 
glm_b_bio_r <- list()
for (i in 1:k) {
  glm_b_bio_r[[i]] <- median(serie_indice_bio_relativa$X50.[serie_indice_bio_relativa$tiempo == i])
}
glm_b_bio_r <- unlist(glm_b_bio_r)

### GAM ---------------------------
gam_b_bio_r <- list()
for (i in 1:k) {
  gam_b_bio_r[[i]] <- median(modelo_gam_r_R2_pred$mu[modelo_gam_r_R2_pred$tiempo == i])
}
gam_b_bio_r <- unlist(gam_b_bio_r)

### Geostatistical ---------------------------
load("./B-Models/relative_biomass/random_geostatistic.RData")
geo_b_bio_r <- list()
for (i in 1:k) {
  geo_b_bio_r[[i]] <- median(m_prj_median[[i]])
}
geo_b_bio_r <- unlist(geo_b_bio_r)

# Save the series of relative biomass ---------------------------
save(gam_b_bio_r, 
     gam_f_bio_r,
     geo_b_bio_r,
     glm_b_bio_r,
     glm_f_bio_r,
     file = "./C-Estimated_series/r_biomass_series.RData"
     )

# CPUE indices ---------------------------
load("./B-Models/CPUE/GLM_GAM_pref.RData")
## Frequentist ---------------------------
### GLM ---------------------------
glm_f_CPUE <- list()
for (i in 1:k) {
  glm_f_CPUE[[i]] <- median(data_pred_glm$respuesta[data_pred_glm$tiempo == i])
}
glm_f_CPUE <- unlist(glm_f_CPUE)

### GAM ---------------------------
gam_f_CPUE <- list()
for (i in 1:k) {
  gam_f_CPUE[[i]] <- median(data_pred_gam$respuesta[data_pred_gam$tiempo == i])
}
gam_f_CPUE <- unlist(gam_f_CPUE)

## Bayesiano ---------------------------
### GLM --------------------------- 
glm_b_CPUE <- list()
for (i in 1:k) {
  glm_b_CPUE[[i]] <- median(serie_indice_bio_relativa$X50.[serie_indice_bio_relativa$tiempo == i])
}
glm_b_CPUE <- unlist(glm_b_CPUE)

### GAM ---------------------------
gam_b_CPUE <- list()
for (i in 1:k) {
  gam_b_CPUE[[i]] <- median(data_gam_r_R2_pred$mu[data_gam_r_R2_pred$tiempo == i])
}
gam_b_CPUE <- unlist(gam_b_CPUE)

### Geostatistical ---------------------------
load("./B-Models/CPUE/pref_geostatistic.RData.RData")
geo_b_CPUE <- list()
for (i in 1:k) {
  geo_b_CPUE[[i]] <- median(m_prj_median[[i]])
}
geo_b_CPUE <- unlist(geo_b_CPUE)

### Marked point process ---------------------------
load("./B-Models/CPUE/marked_point_process.RData")
geo_pref_b_CPUE <- list()
for (i in 1:k) {
  geo_pref_b_CPUE[[i]] <- median(mark$median[mark$tiempo == i])
}
geo_pref_b_CPUE <- unlist(geo_pref_b_CPUE)

# Save series of CPUE indices ---------------------------
save(gam_b_CPUE,
     gam_f_CPUE,
     geo_b_CPUE,
     geo_pref_b_CPUE,
     glm_b_CPUE,
     glm_f_CPUE,
     file = "./C-Series_predichas/CPUE_series.RData"
     )
