## Code to make plots follows.
library(INLA)
library(RColorBrewer)
library(fields)
library(rgdal)
library(raster)
load(paste0("prelim-data", ".RData"))
load(paste0("all-species-output", ".RData"))
load("sighting.RData")
## Species information.
n.species <- length(fit)
species.names <- c("bryde", "cdolp", "bdolp", "orca", "whale", "brydeplus")
## Model information.
n.models <- length(fit[[1]])
model.names <- names(fit[[1]])
## Calculating AICs and determining convergence.
aic.tab <- matrix(0, nrow = n.species, ncol = n.models)
aic.diff.tab <- matrix(0, nrow = n.species, ncol = n.models)
converged.tab <- matrix(FALSE, nrow = n.species, ncol = n.models)
fitted.tab <- matrix(FALSE, nrow = n.species, ncol = n.models)
aic <- function(x) 2*x$objective + 2*length(x$par)
for (i in 1:n.species){
    for (j in 1:n.models){
        if (is.null(fit[[i]][[j]])){
            aic.tab[i, j] <- NA
            converged.tab[i, j] <- NA
        } else {
            fitted.tab[i, j] <- TRUE
            aic.tab[i, j] <- aic(fit[[i]][[j]])
            converged.tab[i , j] <- !any(is.na(rep.summary[[i]][[j]]))
        }
    }
    aic.diff.tab[i, ] <- aic.tab[i, ] - min(aic.tab[i, ], na.rm = TRUE)
}
aic.best.tab <- matrix(FALSE, nrow = n.species, ncol = n.models)
rownames(aic.tab) <- rownames(converged.tab) <- rownames(aic.diff.tab) <- 
    rownames(aic.best.tab) <- rownames(fitted.tab) <- species.names[1:n.species]
colnames(aic.tab) <- colnames(converged.tab) <- colnames(aic.diff.tab) <- 
    colnames(aic.best.tab) <- colnames(fitted.tab) <- model.names
aic.converged.tab <- aic.tab
aic.converged.tab[!converged.tab] <- NA
aic.best.converged.tab <- aic.best.tab

for (i in 1:n.species){
    aic.best.tab[i, which(aic.tab[i, ] == min(aic.tab[i, ], na.rm = TRUE))] <- TRUE
    aic.best.converged.tab[i, which(aic.converged.tab[i, ] ==
                                 min(aic.converged.tab[i, ], na.rm = TRUE))] <- TRUE
}

## Choose a species:
## 1 = "byrde",
## 2 = "cdolp",
## 3 = "bdolp",
## 4 = "orca",
## 5 = "whale"
## 6 = "brydeplus"
## Choose a species.
s <- 1
## Choose a model.
m <- 6
## Grabbing the objects related to this model fit.
fit.use <- fit[[s]][[m]]
d.full.use <- d.full[[s]][[m]]
rep.summary.use <- rep.summary[[s]][[m]]
rand.summary.use <- rand.summary[[s]][[m]]

## Choose month to plot (1 = Aug 2000, 2 = Sep 2000, etc).
i <- 100
## Plotting the spatiotemporal esimates.
proj <- inla.mesh.projector(mesh)
field.proj <- inla.mesh.project(proj, d.full.use[1, , i])

## Set colour scheme with col argument. Set range of z-axis with zmax;
## you'll probably need to change this for plotting estimates for
## individual species.
zmax <- quantile(d.full.use[1, , ], 0.99)
## Choosing a colour scheme for the plots.
cols <- brewer.pal(9, "Blues")
field.proj[field.proj > zmax] <- zmax
## Making a plot.
image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
           zlim = c(0, zmax), main = monthyear.id[i])

## Making a plot of the coastline.
NZ <- readOGR(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
bbox.orig <- data.frame(x = range(obs.xc), y = range(obs.yc))
bbox.small <- bbox.big <- bbox.orig
bbox.small[, 1] <- c(bbox.orig[1, 1] - 0.25*diff(range(bbox.orig[, 1])),
               bbox.orig[2, 1] + 0.25*diff(range(bbox.orig[, 1])))
bbox.small[, 2] <- c(bbox.orig[1, 2] - 0.25*diff(range(bbox.orig[, 2])),
               bbox.orig[2, 2] + 0.25*diff(range(bbox.orig[, 2])))
bbox.big[, 1] <- c(bbox.orig[1, 1] - diff(range(bbox.orig[, 1])),
               bbox.orig[2, 1] + diff(range(bbox.orig[, 1])))
bbox.big[, 2] <- c(bbox.orig[1, 2] - diff(range(bbox.orig[, 2])),
               bbox.orig[2, 2] + diff(range(bbox.orig[, 2])))
s.bbox <- SpatialPoints(bbox.big)
NZ <- crop(NZ, bbox(s.bbox))
## This will create plots of estimated probabilities for every month,
## which can be turned into a .gif. If you uncomment jpeg() and
## dev.off(), it's set up to dump them in /tmp on a Linux
## system. (Feel free to ignore; I can make the gifs easily.)
do.gif <- TRUE
if (do.gif){
    del.files <- list.files("/tmp/hauraki-gifs/")
    file.remove(paste0("/tmp/hauraki-gifs/", del.files))
    dir.create("/tmp/hauraki-gifs/")
    zmax <- quantile(d.full.use[1, , ], 0.99)
    #if (s == 3){
    ##zmax <- 0.05
    #}
    cols <- brewer.pal(9, "Blues")
    for (i in 1:n.months){
        gif.index <- as.character(i)
        zeroes.needed <- nchar(as.character(n.months)) - nchar(gif.index)
        gif.index <- paste0(paste(rep(0, zeroes.needed), collapse = ""), gif.index)
        jpeg(paste0("/tmp/hauraki-gifs/plot", gif.index, ".jpg"))
        field.proj <- inla.mesh.project(proj, d.full.use[1, , i])
        field.proj[field.proj > zmax] <- zmax
        image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
                   zlim = c(0, zmax), main = monthyear.id[i])
        which.counts <- which(y[month.id == i, s] > 0)
        points(obs.xc[which.counts], obs.yc[which.counts], pch = 16, col = "red")
        plot(NZ, col = "grey", add = TRUE)
        cat(i, "of", n.months, "\n")
        dev.off()
    }
}

system("convert -delay 20 -loop 0 /tmp/hauraki-gifs/*.jpg ~/Desktop/heresagif.gif")

## Plotting survey effort.
v <- sighting.df$av.vesselprob
u.x <- sort(unique(obs.xc))
u.y <- sort(unique(obs.yc))
z.v <- matrix(0, nrow = length(u.y), ncol = length(u.x))
for (xi in 1:length(u.x)){
    for (yi in 1:length(u.y)){
        w.p <- which(obs.xc == u.x[xi] & obs.yc == u.y[yi])
        if (length(w.p) == 0){
            z.v[yi, xi] <- 0
        } else {
            z.v[yi, xi] <- v[w.p[1]]
        }
    }
}
z.max <- quantile(z.v, 0.95)
z.v[z.v > z.max] <- z.max
u.x.offset <- diff(u.x[1:2])
u.y.offset <- diff(u.y[1:2])
image.plot(list(x = c(u.x - u.x.offset, u.x[length(u.x)] + u.x.offset),
                y = c(u.y - u.y.offset, u.y[length(u.y)] + u.y.offset),
                z = t(z.v)),
           col = cols, zlim = c(0, z.max), axes = FALSE)
box()
plot(NZ, col = "grey", add = TRUE)

## Plotting changes in occupancy at a specific location over time.

## First, choose a location by selecting a number from the following plot.
field.proj <- inla.mesh.project(proj, d.full[[s]][[m]][1, , i])
image.plot(list(x = proj$x, y = proj$y, z = 1000*field.proj), col = cols,
           zlim = c(0, zmax), main = monthyear.id[i], asp = 1)
#points(obs.xc, obs.yc, pch = ".", cex = 3)
text(mesh$loc[, 1], mesh$loc[, 2], labels = 1:nrow(mesh$loc), cex = 0.75)
plot(NZ, col = "grey", add = TRUE)
## To see the lat/long coordinates of each location, look in mesh$loc.

## Enter location IDs here.
ps <- c(310, 311, 292, 295)
n.plots <- length(ps)
cols <- brewer.pal(6, "Set1")
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
        m <- which(aic.best.tab[s, ])
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
        plot.window(xlim = bbox.small[, 1], ylim = bbox.small[, 2], asp = 1)
        box()
        plot(NZ, col = "grey", add = TRUE)
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

## This is the plot of how sighting probabilities across space are
## affected by changing tempartures (for model int) or
## time-of-year (model int.p). Red locations have increased
## sighting probabilities as temperatures increase, blue areas
## have increased sighting probabilities as temperature decreases
## (for model int). Likewise for model int.p, but the
## increases/decreases are related to the plotted sinusoidal
## function.

## Choose a species. Same codes as above.
s <- 1

## Getting taus.
tau.u.int <- exp(rep.summary[[s]][["int"]][rownames(rep.summary[[s]][["int"]]) == "log_tau_u_int", 1])
tau.u.int.p <- exp(rep.summary[[s]][["int-p"]][rownames(rep.summary[[s]][["int-p"]]) == "log_tau_u_int", 1])
## Side-by-side plots for models int and int-p.
par(mfrow = c(2, 2))
u.int.est <- rand.summary[[s]][["int"]][rownames(rand.summary[[s]][["int"]]) == "u_int_all", 1]/tau.u.int
u.int.est.p <- rand.summary[[s]][["int-p"]][rownames(rand.summary[[s]][["int-p"]]) == "u_int_all", 1]/tau.u.int.p
proj <- inla.mesh.projector(mesh)
field.proj.int <- inla.mesh.project(proj, u.int.est)
field.proj.int.p <- inla.mesh.project(proj, u.int.est.p)
cols <- rev(brewer.pal(11, "RdBu"))
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.int)), col = cols)
plot(NZ, col = "grey", add = TRUE)
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.int.p)), col = cols)
plot(NZ, col = "grey", add = TRUE)
##points(obs.xc, obs.yc, pch = ".")
## Plotting average temperature for each month.
average.month.temp <- numeric(12)
for (i in 1:12){
    average.month.temp[i] <- mean(month.temp[jmonth == i])
}
plot(average.month.temp, type = "l")
## Plotting related sinusoidal function.
gamma <- rep.summary[[s]][["int-p"]][rownames(rep.summary[[s]][["int-p"]]) == "gamma", 1]
xx <- seq(0, 2*pi, length.out = 1000)
yy <- cos(xx - gamma)
plot(xx, yy, type = "l")


## Same for cosfiltered temperature effect.
tau.cf <- exp(rep.summary[[s]][["cf"]][rownames(rep.summary[[s]][["cf"]]) == "log_tau_u_cf", 1])
u.cf.est <- rand.summary[[s]][["cf"]][rownames(rand.summary[[s]][["cf"]]) == "u_cf_all", 1]/tau.cf
proj <- inla.mesh.projector(mesh)
field.proj.cf <- inla.mesh.project(proj, u.cf.est)
cols <- rev(brewer.pal(11, "RdBu"))
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.cf)), col = cols)
plot(NZ, col = "grey", add = TRUE)

## Getting taus.
tau.u.int <- exp(rep.summary[[s]][["cf"]][rownames(rep.summary[[s]][["cf"]]) == "log_tau_u_int", 1])
u.int.est <- rand.summary[[s]][["cf"]][rownames(rand.summary[[s]][["cf"]]) == "u_int_all", 1]/tau.u.int
proj <- inla.mesh.projector(mesh)
field.proj.int <- inla.mesh.project(proj, u.int.est)
field.proj.int.p <- inla.mesh.project(proj, u.int.est.p)
cols <- rev(brewer.pal(11, "RdBu"))
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.int)), col = cols)
plot(NZ, col = "grey", add = TRUE)
