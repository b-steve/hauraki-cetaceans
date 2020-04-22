## Code to make plots follows.
library(INLA)
library(RColorBrewer)
library(fields)
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

## This will create plots of estimated probabilities for every month,
## which can be turned into a .gif. If you uncomment jpeg() and
## dev.off(), it's set up to dump them in /tmp on a Linux
## system. (Feel free to ignore; I can make the gifs easily.)
#jpeg("/tmp/dplot-int%03d.jpg")
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
#dev.off()

## Plotting changes in occupancy at a specific location over time.

## First, choose a location by selecting a number from the following plot.
field.proj <- inla.mesh.project(proj, d.full[[m]][s, , i])
image.plot(list(x = proj$x, y = proj$y, z = 1000*field.proj), col = cols,
           zlim = c(0, zmax), main = monthyear.id[i], asp = 1)
points(obs.xc, obs.yc, pch = ".", cex = 3)
text(mesh$loc[, 1], mesh$loc[, 2], labels = as.character(1:230))
## To see the lat/long coordinates of each location, look in mesh$loc.

## Enter location ID here.
p <- 118
## Run the following to make the plot. Needs tidying (e.g., better x-axis).
p.loc <- mesh$loc[p, ]
plot.new()
plot.window(xlim = c(0, n.months), ylim = c(0, max(d.full[[m]][, p, ])))
for (s in 1:5){
    lines(1:n.months, d.full[[m]][s, p, ], col = s)
}
box()
axis(1)
axis(2)
abline(v = 6 + (0:30)*12, lty = "dotted")
legend("topleft", c("byrde", "cdolp", "bdolp", "orca", "whale"), col = 1:5, lty = rep(1, 5))


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
