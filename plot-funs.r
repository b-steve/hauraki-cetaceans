## A function to calculate AICs from different models for different species.
calc.aics <- function(){
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
    list(aic = aic.tab, diff = aic.diff.tab, converged = converged.tab, fitted = fitted.tab,
         best = aic.best.tab, best.converged = aic.best.converged.tab)
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
    NZ <- crop(NZ, bbox(s.bbox))
    plot(NZ, col = "grey", add = TRUE)
}

## A function to plot a distribution for a given month.
plot.surf <- function(species = 1, model = 1, month = 1, surf = "d", show.obs = FALSE, zlim = NULL, cols = NULL, main = NULL){
    fit <- fit[[species]][[model]]
    if (surf == "d"){
        z.full <- d.full[[species]][[model]]
    }
    rep.summary <- rep.summary[[species]][[model]]
    rand.summary <- rand.summary[[species]][[model]]
    proj <- inla.mesh.projector(mesh)
    field.proj <- inla.mesh.project(proj, z.full[1, , month])
    if (is.null(zlim)){
        zmax <- quantile(z.full[1, , ], 0.99)
        zlim <- c(0, zmax)
    } else {
        zmax <- zlim[2]
    }
    if (is.null(cols)){
        cols <- brewer.pal(9, "Blues")
    }
    field.proj[field.proj > zmax] <- zmax
    if (is.null(main)){
        main <- monthyear.id[month]
    }
    image.plot(list(x = proj$x, y = proj$y, z = field.proj), col = cols,
               zlim = zlim, main = main)
    plot.coast()
    if (show.obs){
        which.counts <- which(y[month.id == month, species] > 0)
        points(obs.xc[which.counts], obs.yc[which.counts], pch = 16, col = "black")
    }
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
        cat(i, "of", n.months, "\n")
        dev.off()
    }
    system(paste0("convert -delay 20 -loop 0 tmp/*.jpg ", file))
    unlink("tmp", recursive = TRUE)
}
