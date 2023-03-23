#*********************************
# SPiCT                          #
# Relative biomass indices       #
# Random sampling                #
#*********************************
# Alba Fuster-Alonso             #
#*********************************

# Packages ---------------------------
library(spict)
library(knitr)
library(formatR)
library(ellipse)
library(ggplot2)
library(corrplot)
library(icesAdvice)

load("./A-Simulacion/simulated_biomass.RData")
# INPUT 1: catch series ---------------------------
## Extract the absolute biomass for each year ---------------------------
biomasa_simulada <- list()
for (i in 1:k) {
  biomasa_simulada[[i]] <- sum(data_variables[[i]]$biomasa)
}
biomasa_simulada <- unlist(biomasa_simulada)
## Set parameters for the production curve ---------------------------
r <- 1.4
K <- sum(biomasa_simulada)/0.5

## Approximation of the catch ---------------------------
capturas_simuladas <- list()
for (i in 1: (k-1)) {
  capturas_simuladas[[i]] <- biomasa_simulada[i] - biomasa_simulada[(i + 1)] +
  (r*biomasa_simulada[i] * (1-(biomasa_simulada[i]/K)))
}

capturas_simuladas[[15]] <-  r*biomasa_simulada[15]*(1-(biomasa_simulada[15]/K)) # Equilibrium for the last year
capturas_simuladas <- unlist(capturas_simuladas)

## Catch ---------------------------
C_sp_i <- data.frame(obsC = capturas_simuladas, timeC = rep(2000:2014))

# INPUT 2: series of relative indices ---------------------------
load("./C-Estimated_series/r_biomass_series.RData") # relative biomass series of indices standardized with different models
load("./C-Estimated_series/CPUE_series.RData") # CPUE series of indices standardized with different models
I_sp_i <- data.frame(ObsI = model, timeI = 2000:2014) # IMPORTANT: change MODEL for a series of indices (we run the script with all the series)

# Joint inputs ---------------------------
inp <- list(
  timeC = C_sp_i$timeC, # Catch
  obsC = C_sp_i$obsC, # Series of catches
  timeI = list(I_sp_i$timeI), # Index of years
  obsI = list(I_sp_i$ObsI) # CPUE indices
)

inp <- check.inp(inp)
inp$dtc
inp$dteuler <- 1/16
inp$ini$logn <- log(2) # Set p (asymmetry parameter Schafer)
inp$phases$logn <- -1

# Fit SPiCT ---------------------------
res <- fit.spict(inp)

# Diagnostic and validation ---------------------------
## Convergence 
res$opt$convergence # 0 convergence  
all(is.finite(res$sd)) # if it is TRUE is okay 

residuals <- calc.osa.resid(res) # residuals of the model
plotspict.diagnostic(residuals) # test for the diagnostic

## Retrospective analysis
r <- fit.spict(inp) 
rep=retro(r, nretroyear=3) # Eliminate 3 years
plotspict.retro(rep) # plot for the retrospective analysis 
m1=mohns_rho(rep, what = c("FFmsy", "BBmsy"));
m1 # Discrepancy between last year and trajectories

## Production curve
plotspict.production(res) # check that the points fall around the curve

## Model robustness (change the initial values)
set.seed(1234)
a=check.ini(inp, ntrials=30 
            ,verbose = FALSE)
a$check.ini$propchng
a$check.ini$resmat # check distance is close to 0 

# Confidence intervals
fit=res # no more than 1 unit
get.par("logBmBmsy", fit, exp=TRUE, CI = 0.9) 
get.par("logFmFmsy", fit, exp=TRUE, CI = 0.9)

# Estimated biomass vs real biomass
get.par("logB", fit, exp=TRUE, CI = 0.9) 
biomass_scpit <- get.par("logB", res, exp = TRUE)[seq(1,length(inp$time),by=16),2]
biomass_scpit <- biomass_scpit[1:15] 
plot(biomass_scpit, col = "red", type = "l")
plot(simulated_biomass, col = "blue", type = "l")
plotspict.biomass(res) 