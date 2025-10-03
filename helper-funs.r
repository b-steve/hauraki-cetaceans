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
    s.bbox <- SpatialPoints(bbox.big)
    NZ <- suppressWarnings(st_crop(NZ, c(xmin = bbox.big[1, 1], xmax = bbox.big[2, 1],
                                         ymin = bbox.big[1, 2], ymax = bbox.big[2, 2])))
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
