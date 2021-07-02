## Loading in packages.
library(TMB)
library(INLA)
library(RColorBrewer)
library(fields)
library(rgdal)
library(raster)
## Loading in data.
load(paste0("prelim-data", ".RData"))
load(paste0("all-species-output", ".RData"))
load("sighting.RData")
NZ <- readOGR(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
source("plot-funs.r")

calc.aics()


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
## Plotting estimated spatially varying effect of temperature for a particular species/model.
plot.surf(species = 6, model = 8, surf = "int")
## Plotting estimated spatial effect.
plot.surf(species = 3, model = 5, surf = "omega")

pdf(file = "~/Desktop/hotspot.pdf")
plot.surf(species = 3, model = 5, month = 96, show.obs = TRUE)
plot.surf(species = 3, model = 5, month = 102, show.obs = TRUE)
plot.surf(species = 3, model = 5, surf = "int")
dev.off()


## Make a distribution gif for a species-model combination.
save.gif(species = 3, model = 5, show.obs = TRUE)
## Plot survey effort.
plot.effort()

## For example, here's an interaction plot for all species.
par(mfrow = c(3, 2))
best.mod <- c(8, 8, 5, 5, 8, 8)
for (i in 1:6){
    plot.surf(species = i, model = best.mod[i], surf = "int")
}

do.occ <- FALSE
if (do.occ){
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
}


do.int <- FALSE
if (do.int){
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
}

pdf("~/Desktop/st-plots")

do.st.add <- TRUE
if (do.st.add){
    ## Plots of the spatial omega field and temporal psi field for the
    ## st-add models.
    s <- 4
    omega.s.est <- rand.summary[[s]][["st-add"]][rownames(rand.summary[[s]][["st-add"]]) == "omega_s_all", 1]
    proj <- inla.mesh.projector(mesh)
    field.proj.omega <- inla.mesh.project(proj, omega.s.est)
    cols <- rev(brewer.pal(11, "RdBu"))
    image.plot(list(x = proj$x, y = proj$y, z = exp(field.proj.omega)), col = cols, main = species.names[s])
    plot(NZ, col = "grey", add = TRUE)

    
    psi.t.est <- rand.summary[[s]][["st-add"]][rownames(rand.summary[[s]][["st-add"]]) == "psi_t_all", 1]
    plot(psi.t.est, type = "l", main = species.names[s])
    abline(v = 12*(0:20), lty = "dotted")
}

