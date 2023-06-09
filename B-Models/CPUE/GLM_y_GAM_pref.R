#*********************************
# GLM and GAM	                   #
# CPUE data                      #
# Preferential sampling          #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Load data ---------------------------
load("./A-Simulation/simulated_data.RData")

# Packages ---------------------------
library(mgcv)
library(INLA)
library(R2BayesX)
library(tidyverse)
library(ggplot2)
library(inlabru)

# Modify data ---------------------------
data_pref_final$tiempo <- as.factor(data_pref_final$tiempo)
data_pref_final$batimetria_2 <- as.vector(I(data_pref_final$batimetria^2))

# Frequentist ---------------------------
## GLM ---------------------------
modelo_f_r <- glm(CPUE ~ batimetria + batimetria_2 + tiempo,
                  data = data_pref_final,
                  family = Gamma(link="log")
                  )

## Prediction ---------------------------
pred_glm <- predict(modelo_f_r,
                    type = "response"
                    )
pred_glm_intervals <- predict(modelo_f_r,
                              type = "response",
                              se.fit = TRUE)$se.fit

intervalo_upr <- pred_glm + (2 * pred_glm_intervals)
intervalo_lwr <- pred_glm - (2 * pred_glm_intervals)

data_pred_glm <- data.frame(respuesta = as.vector(pred_glm),
                            intervalo_lwr = as.vector(intervalo_lwr),
                            intervalo_upr = as.vector(intervalo_upr),
                            tiempo = rep(1:15, each = 50)
                            )

## GAM ---------------------------
modelo_fgam_r <- gam(CPUE ~ s(batimetria) + te(xcoord,ycoord) + tiempo,
                     data = data_pref_final,
                     family = Gamma(link = "log")
                     )

## Prediction ---------------------------
pred_gam <- predict(modelo_fgam_r,
                    type = "response"
                    )
pred_gam_intervals <- predict(modelo_fgam_r, type = "response", se.fit = TRUE)$se.fit

intervalo_upr <- pred_gam + (2 * pred_gam_intervals)
intervalo_lwr <- pred_gam - (2 * pred_gam_intervals)

data_pred_gam <- data.frame(respuesta = as.vector(pred_gam),
                            intervalo_lwr = as.vector(intervalo_lwr),
                            intervalo_upr = as.vector(intervalo_upr),
                            tiempo = rep(1:15, each = 50)
                            )

# Bayesian ---------------------------
## GLM ---------------------------
modelo_bglm_r <- inla(CPUE ~ batimetria + batimetria_2 + tiempo,
                      family = "Gamma",
                      control.predictor = list(compute = TRUE),
                      control.compute = list(
                        dic = TRUE,
                        waic = TRUE,
                        cpo = TRUE,
                        config = TRUE
                      ),
                      data = data_pref_final
                      )

## Prediction ---------------------------
r <- modelo_bglm_r
r.samples <- inla.posterior.sample(1000, r)
psam <- sapply(r.samples, function(x) {
  intercept <- x$latent %>% rownames(.) %>% stringr::str_detect("^\\(Intercept\\)") %>% x$latent[.,]
  beta_y <- x$latent %>% rownames(. ) %>% stringr::str_detect("^tiempo") %>% x$latent[.,]
  beta_2_pred <- x$latent %>% rownames(. ) %>% stringr::str_detect("^batimetria") %>% x$latent[.,]
  beta_2_pred <- beta_2_pred[1]
  beta_3_pred <- x$latent %>% rownames(. ) %>% stringr::str_detect("^batimetria_2") %>% x$latent[.,]
  
  
  predictor1 <- intercept + 
                beta_2_pred*data_pref_final$batimetria[data_pref_final$tiempo == 1] +
                beta_3_pred*data_pref_final$batimetria_2[data_pref_final$tiempo == 1]
  
  pre=list();l=length(beta_y)
  for (i in 1:l){
    pre[[i]]=intercept +
             beta_y[i] + 
             beta_2_pred*data_pref_final$batimetria[data_pref_final$tiempo == (i+1)] +
             beta_3_pred*data_pref_final$batimetria_2[data_pref_final$tiempo == (i+1)]
  }
  
  predictor=predictor1
  
  for (i in 1:l){
    predictor <- c(predictor, pre[[i]])
  }
  
  exp(predictor)
})


q.sam_al_a <- apply(psam, 1, quantile,
                    c(.025, 0.05, 0.5, 0.95, .975), na.rm =TRUE)

serie_indice_bio_relativa <- data.frame(t(q.sam_al_a))
serie_indice_bio_relativa$tiempo <- rep(1:15, each = 50)

## GAM ---------------------------
## R2BayesX ---------------------------
modelo_gam_r_R2 <- bayesx(CPUE ~ sx(batimetria) + sx(xcoord, ycoord, bs = "te") + tiempo,
                          family = "Gamma",
                          method = "MCMC",
                          data = data_pref_final
                          )


### Prediction ---------------------------
modelo_gam_r_R2_pred <- predict(modelo_gam_r_R2,
                                type = "response"
                                )
data_gam_r_R2_pred <- data.frame(media_bio = as.vector(modelo_gam_r_R2_pred),
                                 tiempo = rep(1:15, each = 50)
                                 )

# Save models ---------------------------
save.image("./B-Modelos/CPUE/GLM_GAM_pref.RData")
