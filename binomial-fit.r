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
do.int.p <- TRUE
do.summary <- FALSE

## Loading in the data.
load("sighting.RData")

## Month labels.
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
## Year labels.
years <- 2000:2019
## Creating month.id variable (1 = Aug 2000, 2 = Sep 2000, etc).
month.id <- numeric(nrow(new.df3))
current.id <- 1
monthyear.id <- NULL
jmonth <- NULL
for (i in years){
    for (j in months){
        month.id[new.df3$year == i & new.df3$month == j] <- current.id
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
## Centring ssts and month.temp for numerical stability. Exactly how
## they are centred doesn't affect results.
ssts.centred <- ssts - mean(month.temp)
month.temp.centred <- month.temp - mean(month.temp)
## IMPORTANT: ssts.centred[i] must match
## month.temp.centred[month.id[i]] for all i. For example:
ssts.centred[150000]; month.temp.centred[month.id[150000]]
## Same for ssts and month.temp:
ssts[200000]; month.temp[month.id[200000]]
## Total number of months.
n.months <- max(month.id)


## Making spatial mesh.
load("pixelcoord.RData")
xc.unique <- pixel_coord$Long
yc.unique <- pixel_coord$Lat
pixel.id <- new.df3$pixel_id
obs.xc <- xc.unique[pixel.id]
obs.yc <- yc.unique[pixel.id]
mesh <- inla.mesh.2d(cbind(obs.xc, obs.yc), min.angle = 24, max.edge = c(0.05, 3), cutoff = 0.05)
spde <- inla.spde2.matern(mesh, alpha = 2)
n.meshnodes <- mesh$n

## Compiling the TMB code.
compile("binomial_fit.cpp")
dyn.load(dynlib("binomial_fit"))

## Setting up data for TMB.

## Total number of rows.
n <- nrow(new.df3)
## Number of trips with sightings for each row and each species.
y <- as.matrix(new.df3[, c(8:11, 15)])
## Number of species.
n.species <- ncol(y)
## Number of trips for each row.
n.trials <- new.df3$total.trips

## A design matrix for monthly sighting probability estimates. This
## needs one value per month, starting at August 2000 and ending June
## 2019. The one below only uses temperature.
mat.pred <- cbind(1, month.temp)
## This one is set up for periodic regression.
mat.pred.p <- cbind(1, sin(month.jmonth.rad), cos(month.jmonth.rad))

## Setting up the full design matrix (one row for each observation
## here). Whatever goes in here affects how many animals are estimated
## to be in the gulf in general---a nonspatial effect.
mat <- mat.pred[month.id, ]
## ... And for periodic regression.
mat.p <- mat.pred.p[month.id, ]
## Cell visitation probabilities.
v <- new.df3$av.vesselprob
## Putting it all in a list.
data <- list(n = n, y = y, n_species = n.species, n_trials = n.trials,
             n_betas = ncol(mat), mat = mat,
             ssts = ssts, ssts_centred = ssts.centred,
             jmonth_rad = jmonth.rad, month_jmonth_rad = month.jmonth.rad,
             v = v, month_id = month.id - 1, mesh_id = mesh$idx$loc - 1,
             n_months = n.months, n_meshnodes = n.meshnodes,
             month_temp_centred = month.temp.centred,
             mat_pred = mat.pred,
             spde = spde$param.inla[c("M0","M1","M2")],
             fit_st = 0, fit_int = 0)
## Parameters for TMB.
parameters <- list(betas = matrix(0, nrow = n.species, ncol = ncol(mat)),
                   link_phi = numeric(n.species),
                   log_sigma_u_t = numeric(n.species),
                   log_kappa_u_s = numeric(n.species),
                   log_kappa_u_int = numeric(n.species),
                   log_tau_u_int = numeric(n.species),
                   link_gamma = numeric(n.species),
                   u_st_all = array(0, dim = c(n.species, n.meshnodes, n.months)),
                   u_int_all = matrix(0, nrow = n.species, ncol = n.meshnodes))
## Data for periodic regression.
data.p <- data
data.p$n_betas <- ncol(mat.p)
data.p$mat <- mat.p
data.p$mat_pred <- mat.pred.p
## Parameters for periodic regression.
parameters.p <- parameters
parameters.p$betas <- matrix(0, nrow = n.species, ncol = ncol(mat.p))

if (do.fixed){
    ## Making TMB object for fixed-effects only model. This model has no
    ## spatiotemporal effects. It only allows sighting probabilities to
    ## depend on SST and visitation probabilities. This model should take
    ## less than 1 second to fit.
    obj.fixed <- MakeADFun(data = data, parameters = parameters,
                           map = list(u_st_all = factor(rep(NA, length(parameters$u_st))),
                                      u_int_all = factor(rep(NA, length(parameters$u_int))),
                                      link_phi = factor(rep(NA, length(parameters$link_phi))),
                                      log_sigma_u_t = factor(rep(NA, length(parameters$link_phi))),
                                      log_kappa_u_s = factor(rep(NA, length(parameters$log_kappa_u_s))),
                                      log_kappa_u_int = factor(rep(NA, length(parameters$log_kappa_u_int))),
                                      log_tau_u_int = factor(rep(NA, length(parameters$log_tau_u_int))),
                                      link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                           DLL = "binomial_fit")
    ## Fitting the model.
    fit.fixed <- nlminb(obj.fixed$par, obj.fixed$fn, obj.fixed$gr)
    ## Getting sdreport.
    sdrep.fixed <- sdreport(obj.fixed)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.fixed <- plogis(obj.fixed$report()$d_full_logit)
    ## Saving the fixed-effects model.
    save(fit.fixed, sdrep.fixed, obj.fixed, d.full.fixed, file = "fit-fixed.RData")
} else {
    load("fit-fixed.RData")
}

if (do.fixed.p){
    ## Making TMB object for fixed-effects only model with peridic
    ## regression. This model has no spatiotemporal effects. It only
    ## allows sighting probabilities to depend on month of year and
    ## visitation probabilities. This model should take less than 1
    ## second to fit.
    obj.fixed.p <- MakeADFun(data = data.p, parameters = parameters.p,
                             map = list(u_st_all = factor(rep(NA, length(parameters.p$u_st))),
                                        u_int_all = factor(rep(NA, length(parameters.p$u_int))),
                                        link_phi = factor(rep(NA, length(parameters.p$link_phi))),
                                        log_sigma_u_t = factor(rep(NA, length(parameters.p$link_phi))),
                                        log_kappa_u_s = factor(rep(NA, length(parameters.p$log_kappa_u_s))),
                                        log_kappa_u_int = factor(rep(NA, length(parameters.p$log_kappa_u_int))),
                                        log_tau_u_int = factor(rep(NA, length(parameters.p$log_tau_u_int))),
                                        link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                             DLL = "binomial_fit")
    ## Fitting the model.
    fit.fixed.p <- nlminb(obj.fixed.p$par, obj.fixed.p$fn, obj.fixed.p$gr)
    ## Getting sdreport.
    sdrep.fixed.p <- sdreport(obj.fixed.p)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.fixed.p <- plogis(obj.fixed.p$report()$d_full_logit)
    ## Saving the fixed-effects model.
    save(fit.fixed.p, sdrep.fixed.p, obj.fixed.p, d.full.fixed.p, file = "fit-fixed-p.RData")
} else {
    load("fit-fixed-p.RData")
}

## Making TMB object for spatiotemporal model. This adds a wiggly
## spatial field that varies over time, accounting for spatial and
## temporal correlations in sightings. This model will take somewhere
## between 30 mins and 2 hours to fit, at a guess.
parameters$betas <- matrix(fit.fixed$par, nrow = n.species, ncol = ncol(mat))
data$fit_st <- 1
if (do.st){
    obj.st <- obj.int <- MakeADFun(data = data,
                                   parameters = parameters,
                                   random = "u_st_all",
                                   map = list(u_int_all = factor(rep(NA, length(parameters$u_int))),
                                              log_kappa_u_int = factor(rep(NA, length(parameters$log_kappa_u_int))),
                                              log_tau_u_int = factor(rep(NA, length(parameters$log_tau_u_int))),
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
    save(fit.st, sdrep.st, obj.st, d.full.st, file = "fit-st.RData")
} else {
    load("fit-st.RData")
}

## Making TMB object for spatiotemporal model, similar to above, but
## by using a sinosoidal function of time-of-year instead of
## temperature. This model will take somewhere between 30 mins and 2
## hours to fit, at a guess.
parameters.p$betas <- matrix(fit.fixed.p$par, nrow = n.species, ncol = ncol(mat))
data.p$fit_st <- 1
if (do.st.p){
    obj.st.p <- MakeADFun(data = data.p,
                          parameters = parameters.p,
                          random = "u_st_all",
                          map = list(u_int_all = factor(rep(NA, length(parameters.p$u_int))),
                                     log_kappa_u_int = factor(rep(NA, length(parameters.p$log_kappa_u_int))),
                                     log_tau_u_int = factor(rep(NA, length(parameters.p$log_tau_u_int))),
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
    save(fit.st.p, sdrep.st.p, obj.st.p, d.full.st.p, file = "fit-st-p.RData")
} else {
    load("fit-st-p.RData")
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
parameters$link_phi <- fit.st$par[names(fit.st$par) == "link_phi"]
parameters$log_sigma_u_t <- fit.st$par[names(fit.st$par) == "log_sigma_u_t"]
parameters$log_kappa_u_s <- fit.st$par[names(fit.st$par) == "log_kappa_u_s"]
data$fit_st <- 1
data$fit_int <- 1
if (do.int){
    obj.int <- MakeADFun(data = data,
                         parameters = parameters,
                         random = c("u_st_all", "u_int_all"),
                         map = list(link_gamma = factor(rep(NA, length(parameters$link_gamma)))),
                         inner.control = list(maxit = 50),
                         DLL = "binomial_fit")
    ## Fitting the model.
    fit.int <- nlminb(obj.int$par, obj.int$fn, obj.int$gr)
    ## Getting the sdreport.
    sdrep.int <- sdreport(obj.int)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.int <- plogis(obj.int$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    save(fit.int, sdrep.int, obj.int, d.full.int, file = "fit-int.RData")
} else {
    load("fit-int.RData")
}

## Making TMB object for spatiotemporal model with a
## spatial-temperature interaction, similar to above, but by using a
## sinosoidal function of time-of-year instead of temperature for both
## the fixed and interaction effects. This model will take somewhere
## between 30 mins and 2 hours to fit, at a guess.
parameters.p$betas <- matrix(fit.st.p$par[names(fit.st.p$par) == "betas"], nrow = n.species, ncol = ncol(mat))
parameters.p$link_phi <- fit.st.p$par[names(fit.st.p$par) == "link_phi"]
parameters.p$log_sigma_u_t <- fit.st.p$par[names(fit.st.p$par) == "log_sigma_u_t"]
parameters.p$log_kappa_u_s <- fit.st.p$par[names(fit.st.p$par) == "log_kappa_u_s"]
data.p$fit_st <- 1
data.p$fit_int <- 2
if (do.int.p){
    obj.int.p <- MakeADFun(data = data.p,
                         parameters = parameters.p,
                         random = c("u_st_all", "u_int_all"),
                         map = list(link_gamma = factor(rep(NA, length(parameters.p$link_gamma)))),
                         inner.control = list(maxit = 50),
                         DLL = "binomial_fit")
    ## Fitting the model.
    fit.int.p <- nlminb(obj.int.p$par, obj.int.p$fn, obj.int.p$gr)
    ## Getting the sdreport.
    sdrep.int.p <- sdreport(obj.int.p)
    ## Calculating estimates of sighting probabilities given visitation.
    d.full.int.p <- plogis(obj.int.p$report()$d_full_logit)
    ## Saving the full model with the temperature-interaction field.
    save(fit.int.p, sdrep.int.p, obj.int.p, d.full.int.p, file = "fit-int-p.RData")
} else {
    load("fit-int-p.RData")
}

if (do.summary){
    ## Loading in models, if they've been fitted in a previous R session.
    load("fit-fixed.RData")
    load("fit-st.RData")
    load("fit-int.RData")
    
    ## Comparing models by AIC.
    2*fit.fixed$objective + 2*length(fit.fixed$par) ## Fixed-effects model.
    2*fit.st$objective + 2*length(fit.st$par) ## Spatiotemporal model.
    2*fit.int$objective + 2*length(fit.int$par) ## Spatiotemporal model with space-temperature interaction.

    ## Choose a species.
    s <- 2

    ## Select a model.
    obj <- obj.st
    sdrep <- sdrep.st
    d.full <- d.full.st[s, , ]
    ## Collecting random and reported summaries.
    rand.summary <- summary(sdrep, type = "random")
    rep.summary <- summary(sdrep, type = "report")
    
    ## Plotting the spatiotemporal esimates.
    proj <- inla.mesh.projector(mesh)
    ## Choose month to plot (1 = Aug 2000, 2 = Sep 2000, etc).
    i <- 1
    field.proj <- inla.mesh.project(proj, d.full[, i])
    
    ## Set colour scheme with col argument. Set range of z-axis with zmax;
    ## you'll probably need to change this for plotting estimates for
    ## individual species.
    zmax <- quantile(d.full, 0.99)
    ## Choosing a colour scheme for the plots.
    cols <- brewer.pal(9, "Blues")
    field.proj[field.proj > zmax] <- zmax
    ## Making a plot.
    image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
               zlim = c(0, zmax), main = monthyear.id[i])
    
    ## This will create a plot of estimated probabilities for every
    ## month. It's set up to dump them in /tmp on a Linux system.
    jpeg("/tmp/dplot-int%03d.jpg")
    zmax <- quantile(d.full, 0.99)
    cols <- brewer.pal(9, "Blues")
    for (i in 1:n.months){
        field.proj <- inla.mesh.project(proj, d.full[, i])
        field.proj[field.proj > zmax] <- zmax
        image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
                   zlim = c(0, zmax), main = monthyear.id[i])
        points(obs.xc, obs.yc, pch = ".")
        cat(i, "of", n.months, "\n")
    }
    dev.off()
    
    ## This is the plot of how sighting probabilities across space are
    ## affected by changing tempartures. Red locations have increased
    ## sighting probabilities as temperatures increase, blue areas have
    ## increased sighting probabilities as temperature decreases.
    u.int.est <- rand.summary[rownames(rand.summary) == "u_int", 1]
    proj <- inla.mesh.projector(mesh)
    field.proj <- inla.mesh.project(proj, u.int.est)
    cols <- rev(brewer.pal(11, "RdBu"))
    image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj)), col = cols)
    points(obs.xc, obs.yc, pch = ".")
    
    save.image("fit-everything.RData")
}
