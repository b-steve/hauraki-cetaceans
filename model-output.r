## Loading in packages.
library(TMB)
library(INLA)
library(RColorBrewer)
library(fields)
library(sp)
library(sf)
## Loading in data.
load("prelim-data.RData")
load("all-species-output.RData")
load("sighting.RData")
NZ <- read_sf(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
source("helper-funs.r")
## Setting up sighting counts for the different species categories.
y <- as.matrix(sighting.df[, c(8:11, 15)])
## The final column is the Bryde's plus general whale categories.
y <- cbind(y, y[, 1] + y[, 5])

## The following data frame indicates which effects appear in which
## models.
model.df

## AICs for model selection.
aics <- calc.aics()
## This is a list with a bunch of components:
## - $aic is the AIC value for each species-model pair.
## - $diff is the difference between a model's AIC and the model with
##   the best AIC for that species, and reflects the AIC table in the
##   manuscript.
## - $best has a TRUE entry in the best model (by AIC) for each
##   species.
aics

## A vector of species names.
s.names <- c("Bryde's whale", "Common dolphin", "Bottlenose dolphin",
             "Killer Whale", "Whale", "Bryde's whale +")
## For lots of the codes below, use these values for the species argument:
## 1 = "byrde",
## 2 = "cdolp",
## 3 = "bdolp",
## 4 = "orca",
## 5 = "whale"
## 6 = "brydeplus"

## And choose a model number from model.df

## Plotting estimated relative density for a particular month/species/model.
plot.surf(species = 1, model = 21, month = 123)
## Plotting estimated spatially varying effect of temperature for a
## particular species/model. Only makes sense for models that include
## such an effect, where the 'spatially.varying' column is "SST" in
## model.df.
plot.surf(species = 1, model = 16, surf = "int")
## Plotting estimated spatial effect. This only makes sense for models
## that have a single constant wiggly spatial field over time, where
## the 'spatiotemporal' column is "Space" or "Space+Time" in model.df.
plot.surf(species = 1, model = 14, surf = "omega")

## Plotting estimated temporal effect. This only makes sense for
## models that have a wiggly temporal process over time, where the
## 'spatiotemporal' column is "Time" or "Space+Time" in model.df.
plot.temporal(species = 1, model = 16)

## Make a distribution gif for a species/model combination. The code
## below creates a gif for each species, using the best model for that
## species. Takes a while to run. Note that this function requires the
## command-line tool convert. You can check by running the following
## in R:
## system("convert --version")
## If available, you will see output related to your version of
## convert. The gif is saved to the current R working directory.
for (i in 1:6) save.gif(species = i, model = which(aics$best[i, ]))

## Plot survey effort.
plot.effort()

## For example, here's a plot of the spatially varying effect of
## temperature for all species.
#par(mfrow = c(3, 2))
## Just using the model 16 for these, which isn't actually the best
## model for every species.
#for (i in 1:6){
#    plot.surf(species = i, model = 16, surf = "int")
#}

## Plotting changes in occurrence probability at a specific location
## over time.
    
## Uncomment and run the code below to choose location(s)
## by selecting a number from the following plot. You might need
## to make the graphics window larger to see all the numbers.
#plot.surf(species = 1, model = 1, month = 1)
#text(mesh$loc[, 1], mesh$loc[, 2], labels = 1:nrow(mesh$loc), cex = 0.75)
## To see the lat/long coordinates of each location, look in mesh$loc.

## The locations IDs from above go here.
ps <- c(310, 311, 292, 295)
n.plots <- length(ps)
cols <- brewer.pal(7, "Set1")[-6][order(order(s.names))]
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
par(mar = c(0, 0, 0.25, 0.5), oma = c(6, 6, 0, 0.25), las = 1, xaxs = "i")
all.ds <- array(0, dim = c(n.plots, n.months, 6))
for (i in 1:n.plots){
    p <- ps[i]
    for (s in order(s.names)){
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
    for (s in order(s.names)){
        lines(1:n.months, all.ds[i, , s], col = cols[s])
    }
    box()
    if (i == n.plots){
        axis(1, at = 6 + (0:19)*12, labels = FALSE)
        axis(1, labels = as.character(2000:2019), at = 12*(0:19), tick = FALSE, las = 3)
        mtext("Year", side = 1, line = 4, srt = 90, las = 0)
        mtext("p(s, t)", side = 2, line = 4, las = 0, outer = TRUE)
    }
    axis(2)
    if (i == 1){
        legend("topright", sort(s.names),
               col = cols[order(s.names)], lty = rep(1, 6), bg = "white")
    }
    if (separate.maps){
        plot.new()
        plot.window(xlim = range(proj$x), ylim = range(proj$y), asp = 1)
        box()
        plot.coast()
        pts.cols <- rep(cols[1], 4)
        pts.cols[i] <- cols[3]
        points(mesh$loc[ps, 1], mesh$loc[ps, 2], col = pts.cols, pch = 16, cex = 2)
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

## Plotting a spatial covariance function.
for (s in order(s.names)){
    m <- which(aics$best[s, ])
    kappa <- exp(fit[[s]][[m]]$par["log_kappa_epsilon"])
    dd <- seq(0, 1, length.out = 1000)
    yy <- matern.cov(dd, 1, kappa)
    lat.km <- 110.574
    if (s == order(s.names)[1]){
        plot(lat.km*dd, yy, type = "l", xlab = "Distance (km)", ylab = "Correlation", ylim = c(0, 1), col = cols[s])
    } else {
        lines(lat.km*dd, yy, col = cols[s])
    }
}
legend("topright", legend = sort(s.names), col = cols[order(s.names)], lty = rep(1, 6))

## Plotting a temporal covariance function.
for (s in order(s.names)){
    m <- which(aics$best[s, ])
    link.phi <- fit[[s]][[m]]$par["link_phi_epsilon"]
    phi <- 2*exp(link.phi)/(1 + exp(link.phi)) - 1
    tt <- 0:12
    yy <- phi^tt
    if (s == order(s.names)[1]){
        plot(tt, yy, ylim = c(-0.5, 1), col = cols[s], xlab = "Time (months)", ylab = "Correlation", type = "b")
    } else {
        points(tt, yy, col = cols[s], type = "b")
    }
    print(phi)
}
abline(h = 0, lty = "dotted")
legend("bottomright", legend = sort(s.names), col = cols[order(s.names)], pch = rep(1, 6))
