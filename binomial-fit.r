## ## Reading in command-line stuff for single-species models.
args <- commandArgs(trailingOnly = TRUE)
species <- as.numeric(args[1])
## ## NA will default to all species.
print(species)
#species <- NA

## Loading packages.
library(INLA)
library(TMB)
library(fields)
library(RColorBrewer)

## Set whether or not to fit various models. If not, need an .RData
## file.
do.fixed <- TRUE
do.fixed.p <- TRUE
do.st <- TRUE
do.st.p <- TRUE
do.int <- TRUE
do.int.psi <- TRUE
do.int.p <- TRUE
do.cf <- TRUE

## Loading in the data.
load("sighting.RData")
## Filliing in NaNs (fix when we have full data).
cosfilt.temp[1:6] <- cosfilt.temp[7]
cosfilt.temp[222:227] <- cosfilt.temp[221]

## Month labels.
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
## Year labels.
years <- 2000:2019
## Creating month.id variable (1 = Aug 2000, 2 = Sep 2000, etc).
month.id <- numeric(nrow(sighting.df))
current.id <- 1
monthyear.id <- NULL
jmonth <- NULL
for (i in years){
    for (j in months){
        month.id[sighting.df$year == i & sighting.df$month == j] <- current.id
        monthyear.id <- c(monthyear.id, paste(j, i))
        jmonth <- c(jmonth, which(months == j))
        current.id <- current.id + 1
    }
}

## Starting at August 2000 and ending in June 2019.
monthyear.id <- monthyear.id[-c(1:7, 235:240)]
jmonth <- jmonth[-c(1:7, 235:240)] 
month.id <- month.id - min(month.id) + 1
## Converting monthly jmonth to radians.
month.jmonth.rad <- 2*pi*(jmonth - 0.5)/12
## Getting values for each observation.
jmonth.rad <- month.jmonth.rad[month.id]

## Creating ssts from month.temp.
ssts <- month.temp[month.id]
ssts.cf <- cosfilt.temp[month.id]
## Centring ssts and month.temp for numerical stability. Exactly how
## they are centred doesn't affect results.
ssts.centred <- ssts - mean(month.temp)
ssts.cf.centred <- ssts.cf - mean(cosfilt.temp)
month.temp.centred <- month.temp - mean(month.temp)
cosfilt.temp.centred <- cosfilt.temp - mean(cosfilt.temp)
## IMPORTANT: ssts.centred[i] must match
## month.temp.centred[month.id[i]] for all i. For example:
ssts.centred[150000]; month.temp.centred[month.id[150000]]
ssts.cf.centred[150000]; cosfilt.temp.centred[month.id[150000]]
## Same for ssts and month.temp:
ssts[200000]; month.temp[month.id[200000]]
ssts.cf[200000]; cosfilt.temp[month.id[200000]]
## Total number of months.
n.months <- max(month.id)


## Making spatial mesh.
load("pixelcoord.RData")
xc.unique <- pixel_coord$Long
yc.unique <- pixel_coord$Lat
pixel.id <- sighting.df$pixel_id
obs.xc <- xc.unique[pixel.id]
obs.yc <- yc.unique[pixel.id]
#mesh <- inla.mesh.2d(cbind(obs.xc, obs.yc), min.angle = 24, max.edge = c(0.05, 3), cutoff = 0.05)
mesh <- inla.mesh.2d(cbind(obs.xc, obs.yc), min.angle = 25, max.edge = c(0.05, 2),
                     cutoff = 0.03, offset = 0.1)
#NZ <- readOGR(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
#plot(mesh)
#plot(NZ, col = "grey", add = TRUE)
spde <- inla.spde2.matern(mesh, alpha = 2)
n.meshnodes <- mesh$n

## Compiling the TMB code.
compile("binomial_fit.cpp")
dyn.load(dynlib("binomial_fit"))

## Setting up data for TMB.

## Total number of rows.
n <- nrow(sighting.df)
## Number of trips with sightings for each row and each species.
y <- as.matrix(sighting.df[, c(8:11, 15)])
## Number of total species.
n.all.species <- ncol(y)
## Number of trips for each row.
n.trials <- sighting.df$total.trips

## A design matrix for monthly sighting probability estimates. This
## needs one value per month, starting at August 2000 and ending June
## 2019. The one below only uses temperature.
mat.pred <- cbind(1, month.temp)
## This one is set up for periodic regression.
mat.pred.p <- cbind(1, sin(month.jmonth.rad), cos(month.jmonth.rad))
## This one is set up for the cos-filtered temperature.
mat.pred.cf <- cbind(1, month.temp, cosfilt.temp)

## Setting up the full design matrix (one row for each observation
## here). Whatever goes in here affects how many animals are estimated
## to be in the gulf in general---a nonspatial effect.
mat <- mat.pred[month.id, ]
## ... And for periodic regression.
mat.p <- mat.pred.p[month.id, ]
## ... And for cos-filtered temperature.
mat.cf <- mat.pred.cf[month.id, ]
## Cell visitation probabilities.
v <- sighting.df$av.vesselprob
## Make alternations for substes of species here.
##species <- 1
if (is.na(species)){
    species <- 1:n.all.species
}
y <- y[, species, drop = FALSE]
n.species <- ncol(y)

save.image(file = "prelim-data-smalltri.RData")

## Putting it all in a list.
data <- list(n = n, y = y, n_species = n.species, n_trials = n.trials,
             n_betas = ncol(mat), mat = mat,
             ssts = ssts, ssts_centred = ssts.centred, ssts_cf_centred = ssts.cf.centred,
             jmonth_rad = jmonth.rad, month_jmonth_rad = month.jmonth.rad,
             v = v, month_id = month.id - 1, mesh_id = mesh$idx$loc - 1,
             n_months = n.months, n_meshnodes = n.meshnodes,
             month_temp_centred = month.temp.centred,
             cosfilt_temp_centred = cosfilt.temp.centred,
             mat_pred = mat.pred,
             spde = spde$param.inla[c("M0","M1","M2")],
             fit_psi = 0, fit_omega = 0, fit_epsilon = 0, fit_int = 0, fit_cf = 0)
## Parameters for TMB.
parameters <- list(betas = matrix(0, nrow = n.species, ncol = ncol(mat)),
                   link_phi_psi = numeric(n.species),
                   log_sigma_psi = numeric(n.species),
                   log_kappa_omega = numeric(n.species),
                   log_sigma_omega = numeric(n.species),
                   link_phi_epsilon = numeric(n.species),
                   log_sigma_epsilon = numeric(n.species),
                   log_kappa_epsilon = numeric(n.species),
                   log_kappa_u_int = numeric(n.species),
                   log_tau_u_int = numeric(n.species),
                   log_kappa_u_cf = numeric(n.species),
                   log_tau_u_cf = numeric(n.species),
                   link_gamma = numeric(n.species),
                   psi_t_all = matrix(0, nrow = n.species, ncol = n.months),
                   omega_s_all = matrix(0, nrow = n.species, ncol = n.meshnodes),
                   epsilon_st_all = array(0, dim = c(n.species, n.meshnodes, n.months)),
                   u_int_all = matrix(0, nrow = n.species, ncol = n.meshnodes),
                   u_cf_all = matrix(0, nrow = n.species, ncol = n.meshnodes))
## Data for periodic regression.
data.p <- data
data.p$n_betas <- ncol(mat.p)
data.p$mat <- mat.p
data.p$mat_pred <- mat.pred.p
## Parameters for periodic regression.
parameters.p <- parameters
parameters.p$betas <- matrix(0, nrow = n.species, ncol = ncol(mat.p))
## Data for models with cosfiltered temperatures.
data.cf <- data
data.cf$n_betas <- ncol(mat.cf)
data.cf$mat <- mat.cf
data.cf$mat_pred <- mat.pred.cf
## Parameters for cosfiltered temperature model.
parameters.cf <- parameters
parameters.cf$betas <- matrix(0, nrow = n.species, ncol = ncol(mat.cf))

## Loading test fits.
#load("test.RData")
#load("fit-int.RData")
if (do.fixed){
    ## Making TMB object for fixed-effects only model. This model has no
    ## spatiotemporal effects. It only allows sighting probabilities to
    ## depend on SST and visitation probabilities. This model should take
    ## less than 1 second to fit.
    obj.fixed <- MakeADFun(data = data, parameters = parameters,
                           map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                      omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                      epsilon_st_all = factor(rep(NA, length(parameters$epsilon_st))),
                                      u_int_all = factor(rep(NA, length(parameters$u_int_all))),
                                      u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                      link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                      log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                      log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                      log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                      link_phi_epsilon = factor(rep(NA, length(parameters$link_phi_epsilon))),
                                      log_sigma_epsilon = factor(rep(NA, length(parameters$link_phi_epsilon))),
                                      log_kappa_epsilon = factor(rep(NA, length(parameters$log_kappa_epsilon))),
                                      log_kappa_u_int = factor(rep(NA, length(parameters$log_kappa_u_int))),
                                      log_tau_u_int = factor(rep(NA, length(parameters$log_tau_u_int))),
                                      log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                      log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                      link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                           DLL = "binomial_fit")
    ## Fitting the model.
    fit.fixed <- nlminb(obj.fixed$par, obj.fixed$fn, obj.fixed$gr)
    ## Getting sdreport.
    sdrep.fixed <- sdreport(obj.fixed)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.fixed <- plogis(obj.fixed$report()$d_full_logit)
    ## Saving the fixed-effects model.
    if (n.species == 1){
        save(fit.fixed, sdrep.fixed, d.full.fixed,
             file = paste0("fit-fixed-species-", species, ".RData"))
    } else {
        save(fit.fixed, sdrep.fixed, d.full.fixed, file = "fit-fixed.RData")
    }
    rm(obj.fixed)
}

if (do.fixed.p){
    ## Making TMB object for fixed-effects only model with periodic
    ## regression. This model has no spatiotemporal effects. It only
    ## allows sighting probabilities to depend on month of year and
    ## visitation probabilities. This model should take less than 1
    ## second to fit.
    obj.fixed.p <- MakeADFun(data = data.p, parameters = parameters.p,
                             map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                        omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                        epsilon_st_all = factor(rep(NA, length(parameters$epsilon_st))),
                                        u_int_all = factor(rep(NA, length(parameters$u_int_all))),
                                        u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                        link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                        log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                        log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                        log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                        link_phi_epsilon = factor(rep(NA, length(parameters.p$link_phi_epsilon))),
                                        log_sigma_epsilon = factor(rep(NA, length(parameters.p$link_phi_epsilon))),
                                        log_kappa_epsilon = factor(rep(NA, length(parameters.p$log_kappa_epsilon))),
                                        log_kappa_u_int = factor(rep(NA, length(parameters.p$log_kappa_u_int))),
                                        log_tau_u_int = factor(rep(NA, length(parameters.p$log_tau_u_int))),
                                        log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                        log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                        link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                             DLL = "binomial_fit")
    ## Fitting the model.
    fit.fixed.p <- nlminb(obj.fixed.p$par, obj.fixed.p$fn, obj.fixed.p$gr)
    ## Getting sdreport.
    sdrep.fixed.p <- sdreport(obj.fixed.p)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.fixed.p <- plogis(obj.fixed.p$report()$d_full_logit)
    ## Saving the fixed-effects model.
    if (n.species == 1){
        save(fit.fixed.p, sdrep.fixed.p, d.full.fixed.p,
             file = paste0("fit-fixed-p-species-", species, ".RData"))
    } else {
        save(fit.fixed.p, sdrep.fixed.p, d.full.fixed.p, file = "fit-fixed-p.RData")
    }
    rm(obj.fixed.p)
}

## Making TMB object for spatiotemporal model. This adds a wiggly
## spatial fieljd that varies over time, accounting for spatial and
## temporal correlations in sightings. This model will take somewhere
## between 30 mins and 2 hours to fit, at a guess.
parameters$betas <- matrix(fit.fixed$par, nrow = n.species, ncol = ncol(mat))
data$fit_epsilon <- 1
if (do.st){
    obj.st <- MakeADFun(data = data,
                        parameters = parameters,
                        random = "epsilon_st_all",
                        map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                   omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                   u_int_all = factor(rep(NA, length(parameters$u_int_all))),
                                   u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                   link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                   log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                   log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                   log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                   log_kappa_u_int = factor(rep(NA, length(parameters$log_kappa_u_int))),
                                   log_tau_u_int = factor(rep(NA, length(parameters$log_tau_u_int))),
                                   log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                   log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                   link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                        inner.control = list(maxit = 50),
                        DLL = "binomial_fit")
    ## Fitting the model.
    fit.st <- nlminb(obj.st$par, obj.st$fn, obj.st$gr)
    ## Getting sdreport.
    sdrep.st <- sdreport(obj.st)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.st <- plogis(obj.st$report()$d_full_logit)
    ## Saving the spatiotemporal model.
    if (n.species == 1){
        save(fit.st, sdrep.st, d.full.st,
             file = paste0("fit-st-species-", species, ".RData"))
    } else {
        save(fit.st, sdrep.st, d.full.st, file = "fit-st.RData")
    }
    rm(obj.st)
}

## Making TMB object for spatiotemporal model, similar to above, but
## by using a sinosoidal function of time-of-year instead of
## temperature. This model will take somewhere between 30 mins and 2
## hours to fit, at a guess.
parameters.p$betas <- matrix(fit.fixed.p$par, nrow = n.species, ncol = ncol(mat.p))
data.p$fit_epsilon <- 1
if (do.st.p){
    obj.st.p <- MakeADFun(data = data.p,
                          parameters = parameters.p,
                          random = "epsilon_st_all",
                          map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                     omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                     u_int_all = factor(rep(NA, length(parameters$u_int_all))),
                                     u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                     link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                     log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                     log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                     log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                     log_kappa_u_int = factor(rep(NA, length(parameters.p$log_kappa_u_int))),
                                     log_tau_u_int = factor(rep(NA, length(parameters.p$log_tau_u_int))),
                                     log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                     log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                     link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                          inner.control = list(maxit = 50),
                          DLL = "binomial_fit")
    ## Fitting the model.
    fit.st.p <- nlminb(obj.st.p$par, obj.st.p$fn, obj.st.p$gr)
    ## Getting sdreport.
    sdrep.st.p <- sdreport(obj.st.p)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.st.p <- plogis(obj.st.p$report()$d_full_logit)
    ## Saving the spatiotemporal model.
    if (n.species == 1){
        save(fit.st.p, sdrep.st.p, d.full.st.p,
             file = paste0("fit-st-p-species-", species, ".RData"))
    } else {
        save(fit.st.p, sdrep.st.p, d.full.st.p, file = "fit-st-p.RData")
    }
    rm(obj.st.p)
}

## Making TMB object for spatiotemporal model with a
## spatial-temperature interaction. This model builds on the previous
## one by adding an additional spatial field that estimates animal
## preference for particular locations depending on temperature. So,
## it's not looking to see if there is a change in the number of
## sightings if the temperature changes, rather, it is looking to see
## if animals tend to be found in different places when it's warmer
## than when it's colder. This model should take between 1 hour and 3
## hours to fit, at a guess.
parameters$betas <- matrix(fit.st$par[names(fit.st$par) == "betas"], nrow = n.species, ncol = ncol(mat))
parameters$link_phi_epsilon <- fit.st$par[names(fit.st$par) == "link_phi_epsilon"]
parameters$log_sigma_epsilon <- fit.st$par[names(fit.st$par) == "log_sigma_epsilon"]
parameters$log_kappa_epsilon <- fit.st$par[names(fit.st$par) == "log_kappa_epsilon"]
data$fit_epsilon <- 1
data$fit_int <- 1
if (do.int){
    obj.int <- MakeADFun(data = data,
                         parameters = parameters,
                         random = c("epsilon_st_all", "u_int_all"),
                         map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                    omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                    u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                    link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                    log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                    log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                    log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                    log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                    log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                    link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                         inner.control = list(maxit = 50),
                         DLL = "binomial_fit")
    ## Fitting the model.
    fit.int <- nlminb(obj.int$par, obj.int$fn, obj.int$gr)
    ## Getting the sdreport.
    sdrep.int <- sdreport(obj.int)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.int <- plogis(obj.int$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    if (n.species == 1){
        save(fit.int, sdrep.int, d.full.int,
             file = paste0("fit-int-species-", species, ".RData"))
    } else {
        save(fit.int, sdrep.int, d.full.int, file = "fit-int.RData")
    }
    rm(obj.int)
}

## Same as above, but with the spatiotemporal field separated into three parts.
parameters$betas <- matrix(fit.int$par[names(fit.int$par) == "betas"], nrow = n.species, ncol = ncol(mat))
parameters$link_phi_epsilon <- fit.int$par[names(fit.int$par) == "link_phi_epsilon"]
parameters$log_sigma_epsilon <- fit.int$par[names(fit.int$par) == "log_sigma_epsilon"]
parameters$log_kappa_epsilon <- fit.int$par[names(fit.int$par) == "log_kappa_epsilon"]
parameters$log_tau_u_int <- fit.int$par[names(fit.int$par) == "log_tau_u_int"]
parameters$log_kappa_u_int <- fit.int$par[names(fit.int$par) == "log_kappa_u_int"]
data$fit_psi <- 1
data$fit_omega <- 0
data$fit_epsilon <- 1
data$fit_int <- 1
if (do.int.psi){
    obj.int.psi <- MakeADFun(data = data,
                         parameters = parameters,
                         random = c("psi_t_all", "epsilon_st_all", "u_int_all"),
                         map = list(omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                    u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                    log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                    log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                    log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                    log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf))),
                                    link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                         inner.control = list(maxit = 50),
                         DLL = "binomial_fit")
    ## Fitting the model.
    fit.int.psi <- nlminb(obj.int.psi$par, obj.int.psi$fn, obj.int.psi$gr)
    ## Getting the sdreport.
    sdrep.int.psi <- sdreport(obj.int.psi)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.int.psi <- plogis(obj.int.psi$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    if (n.species == 1){
        save(fit.int.psi, sdrep.int.psi, d.full.int.psi,
             file = paste0("fit-int-sep-species-", species, ".RData"))
    } else {
        save(fit.int.psi, sdrep.int.psi, d.full.int.psi, file = "fit-int-sep.RData")
    }
    rm(obj.int.psi)
}

## Making TMB object for spatiotemporal model with a
## spatial-temperature interaction, similar to above, but by using a
## sinosoidal function of time-of-year instead of temperature for both
## the fixed and interaction effects. This model will take somewhere
## between 30 mins and 2 hours to fit, at a guess.
parameters.p$betas <- matrix(fit.st.p$par[names(fit.st.p$par) == "betas"], nrow = n.species, ncol = ncol(mat.p))
parameters.p$link_phi_epsilon <- fit.st.p$par[names(fit.st.p$par) == "link_phi_epsilon"]
parameters.p$log_sigma_epsilon <- fit.st.p$par[names(fit.st.p$par) == "log_sigma_epsilon"]
parameters.p$log_kappa_epsilon <- fit.st.p$par[names(fit.st.p$par) == "log_kappa_epsilon"]
data.p$fit_epsilon <- 1
data.p$fit_int <- 2
if (do.int.p){
    obj.int.p <- MakeADFun(data = data.p,
                           parameters = parameters.p,
                           random = c("epsilon_st_all", "u_int_all"),
                           map = list(psi_t_all = factor(rep(NA, length(parameters$psi_t_all))),
                                      omega_s_all = factor(rep(NA, length(parameters$omega_s_all))),
                                      u_cf_all = factor(rep(NA, length(parameters$u_cf_all))),
                                      link_phi_psi = factor(rep(NA, length(parameters$link_phi_psi))),
                                      log_sigma_psi = factor(rep(NA, length(parameters$log_sigma_psi))),
                                      log_kappa_omega = factor(rep(NA, length(parameters$log_kappa_omega))),
                                      log_sigma_omega = factor(rep(NA, length(parameters$log_sigma_omega))),
                                      log_kappa_u_cf = factor(rep(NA, length(parameters$log_kappa_u_cf))),
                                      log_tau_u_cf = factor(rep(NA, length(parameters$log_tau_u_cf)))),
                           inner.control = list(maxit = 50),
                           DLL = "binomial_fit")
    ## Fitting the model.
    fit.int.p <- nlminb(obj.int.p$par, obj.int.p$fn, obj.int.p$gr)
    ## Getting the sdreport.
    sdrep.int.p <- sdreport(obj.int.p)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.int.p <- plogis(obj.int.p$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    if (n.species == 1){
        save(fit.int.p, sdrep.int.p, d.full.int.p,
             file = paste0("fit-int-p-species-", species, ".RData"))
    } else {
        save(fit.int.p, sdrep.int.p, d.full.int.p, file = "fit-int-p.RData")
    }
    rm(obj.int.p)
}

## Making TMB object for spatiotemporal model with two types of
## spatial-temperature interactions. This model builds on fit.int by
## adding an additional spatial field that estimates animal preference
## for particular locations depending on gradual changes
## temperature. This allows us to see how distributions are changing
## due to long-term temperature variation.
parameters.cf$betas[, 1:2] <- matrix(fit.int$par[names(fit.int$par) == "betas"],
                                     nrow = n.species, ncol = ncol(mat))
parameters.cf$link_phi_epsilon <- fit.int$par[names(fit.int$par) == "link_phi_epsilon"]
parameters.cf$log_sigma_epsilon <- fit.int$par[names(fit.int$par) == "log_sigma_epsilon"]
parameters.cf$log_kappa_epsilon <- fit.int$par[names(fit.int$par) == "log_kappa_epsilon"]
parameters.cf$log_kappa_u_int <- fit.int$par[names(fit.int$par) == "log_kappa_u_int"]
parameters.cf$log_tau_u_int <- fit.int$par[names(fit.int$par) == "log_tau_u_int"]
data.cf$fit_epsilon <- 1
data.cf$fit_int <- 1
data.cf$fit_cf <- 1
if (do.cf){
    obj.cf <- MakeADFun(data = data.cf,
                        parameters <- parameters.cf,
                        random = c("epsilon_st_all", "u_int_all", "u_cf_all"),
                        map = list(psi_t_all = factor(rep(NA, length(parameters.cf$psi_t_all))),
                                   omega_s_all = factor(rep(NA, length(parameters.cf$omega_s_all))),
                                   link_phi_psi = factor(rep(NA, length(parameters.cf$link_phi_psi))),
                                   log_sigma_psi = factor(rep(NA, length(parameters.cf$log_sigma_psi))),
                                   log_kappa_omega = factor(rep(NA, length(parameters.cf$log_kappa_omega))),
                                   log_sigma_omega = factor(rep(NA, length(parameters.cf$log_sigma_omega))),
                                   link_gamma = factor(rep(NA, length(parameters.cf$link_gamma)))),
                        inner.control = list(maxit = 50),
                        DLL = "binomial_fit")
    ## Fitting the model.
    fit.cf <- nlminb(obj.cf$par, obj.cf$fn, obj.cf$gr)
    ## Getting the sdreport.
    sdrep.cf <- sdreport(obj.cf)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.cf <- plogis(obj.cf$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    if (n.species == 1){
        save(fit.cf, sdrep.cf, d.full.cf, 
             file = paste0("fit-cf-species-", species, ".RData"))
    } else {
        save(fit.cf, sdrep.cf, d.full.cf, file = "fit-cf.RData")
    }
    rm(obj.cf)
}
