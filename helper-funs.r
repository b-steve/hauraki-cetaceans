## Loading in data.
load("prelim-data.RData")
load("all-species-output.RData")
load("sighting.RData")
NZ <- read_sf(dsn = "./kx-nz-seacoast-poly-SHP", layer = "nz-seacoast-poly")
study.site <- read_sf(dsn = "kx-nz-seacoast-poly-SHP/erase-square", layer = "erase_square2")

## A function to calculate AICs from different models for different species.
calc.aics <- function(){
    ## Species information.
    n.species <- length(fit)
    species.names <- c("bryde", "cdolp", "bdolp", "orca", "whale", "brydeplus")
    ## Model information.
    n.models <- length(fit[[1]])
    model.names <- names(fit[[1]])
    ## Calculating AICs.
    aic.tab <- matrix(0, nrow = n.species, ncol = n.models)
    aic.diff.tab <- matrix(0, nrow = n.species, ncol = n.models)
    fitted.tab <- matrix(FALSE, nrow = n.species, ncol = n.models)
    aic <- function(x) 2*x$objective + 2*length(x$par)
    for (i in 1:n.species){
        for (j in 1:n.models){
            if (is.null(fit[[i]][[j]])){
                aic.tab[i, j] <- NA
            } else {
                fitted.tab[i, j] <- TRUE
                aic.tab[i, j] <- aic(fit[[i]][[j]])
            }
        }
        aic.diff.tab[i, ] <- aic.tab[i, ] - min(aic.tab[i, ], na.rm = TRUE)
    }
    aic.best.tab <- matrix(FALSE, nrow = n.species, ncol = n.models)
    rownames(aic.tab) <- rownames(aic.diff.tab) <- rownames(aic.best.tab) <-
        rownames(fitted.tab) <- species.names[1:n.species]
    colnames(aic.tab) <- colnames(aic.diff.tab) <- colnames(aic.best.tab) <-
        colnames(fitted.tab) <- model.names    
    for (i in 1:n.species){
        aic.best.tab[i, which(aic.tab[i, ] == min(aic.tab[i, ], na.rm = TRUE))] <- TRUE
    }
    list(aic = aic.tab, diff = aic.diff.tab, best = aic.best.tab)
}

plot.coast <- function(){
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
    NZ <- suppressWarnings(st_crop(NZ, c(xmin = bbox.big[1, 1], xmax = bbox.big[2, 1],
                                         ymin = bbox.big[1, 2], ymax = bbox.big[2, 2])))
    study.site <- suppressWarnings(st_crop(study.site, c(xmin = bbox.big[1, 1], xmax = bbox.big[2, 1],
                                                         ymin = bbox.big[1, 2], ymax = bbox.big[2, 2])))
    suppressWarnings(plot(study.site, col = "white", border = "white", add = TRUE))
    suppressWarnings(plot(NZ, col = "grey", add = TRUE))
}

## A function to plot a distribution for a given month.
plot.surf <- function(species = 1, model = 1, month = 1, surf = "d", show.obs = FALSE, zlim = NULL, cols = NULL, main = NULL){
    fit <- fit[[species]][[model]]
    if (surf == "d"){
        z.full <- d.full[[species]][[model]]
        z <- z.full[1, , month]
        if (is.null(zlim)){
            zlim <- c(0, quantile(z.full[1, , ], 0.99))
        }
        if (is.null(cols)){
            cols <- brewer.pal(9, "Blues")
        }
        if (is.null(main)){
            main <- monthyear.id[month]
        }
    } else if (surf == "omega"){
        z <- rand.summary[[species]][[model]][rownames(rand.summary[[species]][[model]]) == "omega_s_all", 1]
        if (is.null(zlim)){
            zlim <- c(0, quantile(z, 0.99))
        }
        if (is.null(cols)){
            cols <- brewer.pal(9, "Blues")
        }        
    } else if (surf == "int"){
        tau.u.int <- exp(rep.summary[[species]][[model]][rownames(rep.summary[[species]][[model]]) == "log_tau_u_int", 1])
        z <- exp(rand.summary[[species]][[model]][rownames(rand.summary[[species]][[model]]) == "u_int_all", 1]/tau.u.int)
        if (is.null(cols)){
            cols <- rev(brewer.pal(11, "RdBu"))
        }
        if (is.null(zlim)){
            zlim <- range(z)
            max.zlim.diff <- max(abs(c(1 - zlim[1], zlim[2] - 1)))
            zlim <- 1 + c(-1, 1)*max.zlim.diff
        }
    } else if (surf == "int-cf"){
        tau.cf <- exp(rep.summary[[species]][[model]][rownames(rep.summary[[species]][[model]]) == "log_tau_u_cf", 1])
        z <- exp(rand.summary[[species]][[model]][rownames(rand.summary[[species]][[model]]) == "u_cf_all", 1]/tau.cf)
        if (is.null(cols)){
            cols <- rev(brewer.pal(11, "RdBu"))
        }
        if (is.null(zlim)){
            zlim <- range(z)
        }
    }
    zmin <- zlim[1]
    zmax <- zlim[2]
    proj <- inla.mesh.projector(mesh)
    field.proj <- inla.mesh.project(proj, z)
    field.proj[field.proj > zmax] <- zmax
    field.proj[field.proj < zmin] <- zmin
    image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
               zlim = zlim, main = main)
    plot.coast()
    if (surf == "d" & show.obs){
        which.counts <- which(y[month.id == month, species] > 0)
        points(obs.xc[which.counts], obs.yc[which.counts], pch = 16, col = "black")
    }
}

plot.temporal <- function(species = 1, model = 1){
    psi.t.est <- rand.summary[[species]][[model]][rownames(rand.summary[[species]][[model]]) == "psi_t_all", 1]
    plot(psi.t.est, type = "l")
    abline(v = 12*(0:20), lty = "dotted")
}

save.gif <- function(species = 1, model = 1, file = paste0("species-", species, "-mod-", model, ".gif"),
                     show.obs = FALSE, zlim = NULL, cols = NULL, main = NULL){
    dir.create("tmp")
    for (i in 1:n.months){
        gif.index <- as.character(i)
        zeroes.needed <- nchar(as.character(n.months)) - nchar(gif.index)
        gif.index <- paste0(paste(rep(0, zeroes.needed), collapse = ""), gif.index)
        jpeg(paste0("tmp/plot", gif.index, ".jpg"))
        plot.surf(species = species, model = model, month = i, show.obs = show.obs,
                  zlim = zlim, cols = cols, main = main)       
        cat("Creating frame", i, "of", n.months, "\n")
        dev.off()
    }
    system(paste0("convert -delay 20 -loop 0 tmp/*.jpg ", file))
    unlink("tmp", recursive = TRUE)
}

plot.effort <- function(cols = NULL){
    if (is.null(cols)){
        cols <- brewer.pal(9, "Blues")
    }
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
    plot.coast()
}

## Matern covariance function with nu = 1.
matern.cov <- function(d, sigma, kappa){
    nu <- 1
    out <- (sigma^2/(2^(nu - 1) * gamma(nu))) * ((kappa * abs(d))^nu) * 
        besselK(kappa * abs(d), nu)
    out[d == 0] <- sigma^2
    out
}

loc.id.selector <- function(){
    plot.surf(species = 1, model = 1, month = 1)
    text(mesh$loc[, 1], mesh$loc[, 2], labels = 1:nrow(mesh$loc), cex = 0.75)
}

plot.temporal.at.locs <- function(loc.ids, models, separate.maps = TRUE, same.scale = FALSE){
    ## A vector of species names.
    s.names <- c("Bryde's whale", "Common dolphin", "Bottlenose dolphin",
                 "Killer Whale", "Whale", "Bryde's whale +")
    n.plots <- length(loc.ids)
    cols <- brewer.pal(7, "Set1")[-6][order(order(s.names))]
    ## Mesh projection.
    proj <- inla.mesh.projector(mesh)
    if (separate.maps){
        mat.layout <- matrix(c(rep(1:(2*n.plots), times = rep(c(3, 1), n.plots))),
                             ncol = 4, byrow = TRUE)
    } else {
        mat.layout <- cbind(matrix(rep(c(1, 3:(n.plots + 1)), each = 3),
                                   ncol = 3, byrow = TRUE),
                            c(2, rep(0, n.plots - 1)))
    }
    opar <- par()
    layout(mat.layout, heights = rep(1, 4), widths = rep(1, 4))
    par(mar = c(0, 0, 0.25, 0.5), oma = c(6, 6, 0, 0.25), las = 1, xaxs = "i")
    all.ds <- array(0, dim = c(n.plots, n.months, 6))
    for (i in 1:n.plots){
        p <- loc.ids[i]
        for (s in order(s.names)){
            m <- models[s]
            all.ds[i, , s] <- d.full[[s]][[m]][1, p, ]
        }
    }
    for (i in 1:n.plots){
        p <- loc.ids[i]
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
            points(mesh$loc[loc.ids, 1], mesh$loc[loc.ids, 2], col = pts.cols, pch = 16, cex = 2)
        } else {
            if (i == 1){
                plot.new()
                plot.window(xlim = bbox.small[, 1], ylim = bbox.small[, 2], asp = 1)
                box()
                plot(NZ, col = "grey", add = TRUE)
                text(mesh$loc[loc.ids, 1], mesh$loc[loc.ids, 2], labels = as.character(1:n.plots))
            }
        }
    }
    suppressWarnings(par(opar))
}

plot.spatial.cov <- function(models){
    s.names <- c("Bryde's whale", "Common dolphin", "Bottlenose dolphin",
                 "Killer Whale", "Whale", "Bryde's whale +")
    cols <- brewer.pal(7, "Set1")[-6][order(order(s.names))]
    for (s in order(s.names)){
        m <- models[s]
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
}

plot.temporal.cov <- function(models){
    s.names <- c("Bryde's whale", "Common dolphin", "Bottlenose dolphin",
                 "Killer Whale", "Whale", "Bryde's whale +")
    cols <- brewer.pal(7, "Set1")[-6][order(order(s.names))]
    for (s in order(s.names)){
        m <- models[s]
        link.phi <- fit[[s]][[m]]$par["link_phi_epsilon"]
        phi <- 2*exp(link.phi)/(1 + exp(link.phi)) - 1
        tt <- 0:12
        yy <- phi^tt
        if (s == order(s.names)[1]){
            plot(tt, yy, ylim = c(-0.5, 1), col = cols[s], xlab = "Time (months)", ylab = "Correlation", type = "b")
        } else {
            points(tt, yy, col = cols[s], type = "b")
        }
    }
    abline(h = 0, lty = "dotted")
    legend("bottomright", legend = sort(s.names), col = cols[order(s.names)], pch = rep(1, 6))
}
