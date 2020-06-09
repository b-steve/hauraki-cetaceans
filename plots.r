## Code to make plots follows.
library(INLA)
library(RColorBrewer)
library(fields)
library(rgdal)
library(raster)
load("all-output.RData")
## Comparing models by AIC.
sapply(fit, function(x) 2*x$objective + 2*length(x$par))

## Choose a model. Best to keep this at 5, because it's the best model.
m <- 5
## Choose a species:
## 1 = "byrde",
## 2 = "cdolp",
## 3 = "bdolp",
## 4 = "orca",
## 5 = "whale"
s <- 2
## Choose month to plot (1 = Aug 2000, 2 = Sep 2000, etc).
i <- 100
## Plotting the spatiotemporal esimates.
proj <- inla.mesh.projector(mesh)
field.proj <- inla.mesh.project(proj, d.full[[m]][s, , i])

## Set colour scheme with col argument. Set range of z-axis with zmax;
## you'll probably need to change this for plotting estimates for
## individual species.
zmax <- quantile(d.full[[m]][s, , ], 0.99)
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
                                        #jpeg("/tmp/dplot-int%03d.jpg")
do.gif <- FALSE
if (do.gif){
    zmax <- quantile(d.full[[m]][s, , ], 0.99)
    cols <- brewer.pal(9, "Blues")
    for (i in 1:n.months){
        field.proj <- inla.mesh.project(proj, d.full[[m]][s, , i])
        field.proj[field.proj > zmax] <- zmax
        image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
                   zlim = c(0, zmax), main = monthyear.id[i])
        points(obs.xc, obs.yc, pch = ".")
        cat(i, "of", n.months, "\n")
    }
}

## Plotting changes in occupancy at a specific location over time.

## First, choose a location by selecting a number from the following plot.
field.proj <- inla.mesh.project(proj, d.full[[m]][s, , i])
image.plot(list(x = proj$x, y = proj$y, z = 1000*field.proj), col = cols,
           zlim = c(0, zmax), main = monthyear.id[i], asp = 1)
points(obs.xc, obs.yc, pch = ".", cex = 3)
text(mesh$loc[, 1], mesh$loc[, 2], labels = as.character(1:230))
## To see the lat/long coordinates of each location, look in mesh$loc.

## Enter location IDs here.
ps <- c(118, 124, 117, 126)
n.plots <- length(ps)
cols <- brewer.pal(5, "Set1")
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
all.ds <- array(0, dim = c(n.plots, n.months, 5))
for (i in 1:n.plots){
    p <- ps[i]
    for (s in 1:5){
        all.ds[i, , s] <- d.full[[m]][s, p, ]
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
    for (s in 1:5){
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
        legend("topright", c("Bryde's whale", "Common dolphin", "Bottlenose dolphin", "Orca", "Whale"),
               col = cols, lty = rep(1, 5), bg = "white")
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
s <- 2
## Getting taus.
tau.u.int <- exp(rep.summary[["int"]][rownames(rep.summary[["int"]]) == "log_tau_u_int", 1])[s]
tau.u.int.p <- exp(rep.summary[["int.p"]][rownames(rep.summary[["int.p"]]) == "log_tau_u_int", 1])[s]
## Side-by-side plots for models int and int.p.
par(mfrow = c(2, 2))
u.int.est <- matrix(rand.summary[["int"]][rownames(rand.summary[["int"]]) == "u_int_all", 1],
                    nrow = 5)[s, ]/tau.u.int
u.int.est.p <- matrix(rand.summary[["int.p"]][rownames(rand.summary[["int.p"]]) == "u_int_all", 1],
                      nrow = 5)[s, ]/tau.u.int.p
proj <- inla.mesh.projector(mesh)
field.proj.int <- inla.mesh.project(proj, u.int.est)
field.proj.int.p <- inla.mesh.project(proj, u.int.est.p)
cols <- rev(brewer.pal(11, "RdBu"))
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.int)), col = cols)
image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.int.p)), col = cols)
                                        #points(obs.xc, obs.yc, pch = ".")
## Plotting average temperature for each month.
average.month.temp <- numeric(12)
for (i in 1:12){
    average.month.temp[i] <- mean(month.temp[jmonth == i])
}
plot(average.month.temp, type = "l")
## Plotting related sinusoidal function.
gamma <- rep.summary[["int.p"]][rownames(rep.summary[["int.p"]]) == "gamma", 1][s]
xx <- seq(0, 2*pi, length.out = 1000)
yy <- cos(xx - gamma)
plot(xx, yy, type = "l")
