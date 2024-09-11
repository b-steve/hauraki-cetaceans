## Loading in packages.
library(TMB)
library(INLA)
library(RColorBrewer)
library(fields)
library(sf)
## Loading in data.
load(paste0("prelim-data", ".RData"))
load(paste0("all-species-output", ".RData"))
load("sighting.RData")
NZ <- read_sf(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
source("helper-funs.r")
## Setting up sighting counts for the different species categories.
y <- as.matrix(sighting.df[, c(8:11, 15)])
## The final column is the Bryde's plus general whale categories.
y <- cbind(y, y[, 1] + y[, 5])

## Note that in a lot of the code below, objects related to this
## particular data set are assumed to be sitting in the global
## environment (e.g., an object 'y' for the sighting counts, an object
## 'fits' provided by the process-output.r code, an object 'mesh' for
## the spatial mesh, and so on). It's not the best programming
## practice (sorry!) but the code was written for this particular data
## set and there's too much specific going on to make the functions
## work in general without putting in a lot of effort. For a new data
## set, I recommend taking a look at the code inside calc.aics() and
## plot.surf(), then adjusting it for your purposes.

## AICs for model selection.
aics <- calc.aics()
## This is a list with a bunch of components:
## - $aic is the AIC value for each species-model pair.
## - $diff is the difference between a model's AIC and the model with
##   the best AIC for that species.
## - $converged indicates whether or not the optimisation algorithm
##   converged.
## - $fitted indicates whether or not we actually tried to fit the
##   model (there's no point fitting some of the complicated models to
##   species with low numbers of sightings).
## - $best has a TRUE entry in the best model (by AIC) for each
##   species.
## - $best.converged can be ignored.
aics

s.names <- c("Bryde's", "Common dolphin", "Bottlenose dolphin", "Orca", "Whale", "Bryde's plus whale")
## Choose a species:
## 1 = "byrde",
## 2 = "cdolp",
## 3 = "bdolp",
## 4 = "orca",
## 5 = "whale"
## 6 = "brydeplus"
## Choose a model via names(fit[[1]])
## Plotting estimated relative density for a particular month/species/model.
plot.surf(species = 6, model = 8, month = 100, show.obs = TRUE)
## Plotting estimated spatially varying effect of temperature for a
## particular species/model. Only makes sense for models that include
## such an effect.
plot.surf(species = 1, model = 16, surf = "int")
## Plotting estimated spatial effect. This only makes sense for model
## that have a single constant wiggly spatial field over time.
plot.surf(species = 3, model = 13, surf = "omega")

## Make a distribution gif for a species/model combination. Here I'm
## making a gif for each species, using the best model for that
## species. Takes a while to run.
for (i in 1:6) save.gif(species = i, model = which(aics$best[i, ]), show.obs = TRUE)
## Plot survey effort.
plot.effort()

## For example, here's a plot of the spatially varying effect of
## temperature for all species.
par(mfrow = c(3, 2))
## Just using the model 16 for these, which isn't actually the best
## model for every species.
for (i in 1:6){
    plot.surf(species = i, model = 16, surf = "int")
}

## Plotting changes in occurrence probability at a specific location
## over time.
do.occ <- FALSE
if (do.occ){
    
    ## Uncomment and run the code below to choose location(s)
    ## by selecting a number from the following plot. You might need
    ## to make the graphics window larger to see all the numbers.
    #plot.surf(species = 1, model = 1, month = 1)
    #text(mesh$loc[, 1], mesh$loc[, 2], labels = 1:nrow(mesh$loc), cex = 0.75)
    ## To see the lat/long coordinates of each location, look in mesh$loc.
    
    ## The locations IDs from above go here.
    ps <- c(310, 311, 292, 295)
    n.plots <- length(ps)
    cols <- brewer.pal(6, "Set1")
    ## Mesh projection.
    proj <- inla.mesh.projector(mesh)
    ## Set to TRUE for same scale.
    same.scale <- FALSE
    ## Sorting out layout.
    separate.maps <- TRUE
    if (separate.maps){
        mat.layout <- matrix(c(rep(1:(2*n.plots), times = rep(c(3, 1), n.plots))),
                             ncol = 4, byrow = TRUE)
    } else {
        mat.layout <- cbind(matrix(rep(c(1, 3:(n.plots + 1)), each = 3),
                                   ncol = 3, byrow = TRUE),
                            c(2, rep(0, n.plots - 1)))
    }
    layout(mat.layout, heights = rep(1, 4), widths = rep(1, 4))
    par(mar = c(0, 0, 0.25, 0.5), oma = c(3, 4, 0, 0.25), las = 1, xaxs = "i")
    all.ds <- array(0, dim = c(n.plots, n.months, 6))
    for (i in 1:n.plots){
        p <- ps[i]
        for (s in 1:6){
            ## Note: probably need to set a species-specific model here.
            m <- which(aics$best[s, ])
            all.ds[i, , s] <- d.full[[s]][[m]][1, p, ]
        }
    }
    for (i in 1:n.plots){
        p <- ps[i]
        p.loc <- mesh$loc[p, ]
        plot.new()
        if (same.scale){
            ylim <- c(0, max(all.ds))
        } else {
            ylim <- c(0, max(all.ds[i, , ]))
        }
        plot.window(xlim = c(1, n.months),
                    ylim = ylim)
        abline(v = 6 + (0:30)*12, col = "lightgrey")
        for (s in 1:6){
            lines(1:n.months, all.ds[i, , s], col = cols[s])
        }
        box()
        if (i == n.plots){
            axis(1, at = 6 + (0:19)*12, labels = FALSE)
            axis(1, labels = as.character(2000:2019), at = 12*(0:19), tick = FALSE)
            title(xlab = "Year")
        }
        axis(2)
        if (i == 1){
            legend("topright", c("Bryde's whale", "Common dolphin", "Bottlenose dolphin", "Orca", "Whale", "Bryde's +"),
                   col = cols, lty = rep(1, 6), bg = "white")
        }
        if (separate.maps){
            plot.new()
            plot.window(xlim = range(proj$x), ylim = range(proj$y), asp = 1)
            box()
            plot.coast()
            pts.cols <- rep(cols[2], 4)
            pts.cols[i] <- cols[1]
            points(mesh$loc[ps, 1], mesh$loc[ps, 2], col = pts.cols, pch = 16)
        } else {
            if (i == 1){
                plot.new()
                plot.window(xlim = bbox.small[, 1], ylim = bbox.small[, 2], asp = 1)
                box()
                plot(NZ, col = "grey", add = TRUE)
                text(mesh$loc[ps, 1], mesh$loc[ps, 2], labels = as.character(1:n.plots))
            }
        }
    }
}

do.st.add <- TRUE
if (do.st.add){
    ## Plots of the spatial omega field and temporal psi field for the
    ## st-add models.
    s <- 1
    m <- 11
    plot.surf(species = s, model = m, surf = "omega")

    ## Plotting temporal process.
    psi.t.est <- rand.summary[[s]][[m]][rownames(rand.summary[[s]][[m]]) == "psi_t_all", 1]
    plot(psi.t.est, type = "l")
    abline(v = 12*(0:20), lty = "dotted")
}

png("epsilon-space.png")
cols <- brewer.pal(6, "Set1")
for (s in 1:6){
    m <- which(aics$best[s, ])
    kappa <- exp(fit[[s]][[m]]$par["log_kappa_epsilon"])
    dd <- seq(0, 1, length.out = 1000)
    yy <- matern.cov(dd, 1, kappa)
    lat.km <- 110.574
    if (s == 1){
        plot(lat.km*dd, yy, type = "l", xlab = "Distance (km)", ylab = "Correlation", ylim = c(0, 1), col = cols[s])
    } else {
        lines(lat.km*dd, yy, col = cols[s])
    }
}
legend("topright", legend = s.names, col = cols, lty = rep(1, 6))
dev.off()

png("epsilon-time.png")
for (s in 1:6){
    m <- which(aics$best[s, ])
    link.phi <- fit[[s]][[m]]$par["link_phi_epsilon"]
    phi <- 2*exp(link.phi)/(1 + exp(link.phi)) - 1
    tt <- 0:12
    yy <- phi^tt
    if (s == 1){
        plot(tt, yy, ylim = c(-0.5, 1), col = cols[s], xlab = "Time (months)", ylab = "Correlation")
    } else {
        points(tt, yy, col = cols[s])
    }
    print(phi)
}
abline(h = 0, lty = "dotted")
legend("bottomright", legend = s.names, col = cols, pch = rep(1, 6))
dev.off()
