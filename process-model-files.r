library(TMB)
## Number of species.
n.species <- 6
## Model names from model fitting code.
model.ext <- c("fixed", "fixed-p", "sv", "sv-p", "fixed-s", "fixed-s-p", "fixed-s-int", "fixed-s-int-p",
               "fixed-t", "fixed-t-p", "fixed-t-int", "fixed-t-int-p", "st-add", "st-add-p", "st-add-int",
               "st-add-int-p", "st-nofixed", "st", "st-p", "int-nofixed", "int", "int-p-sst")

## R object indicator for each of the models.
model.names <- chartr("-", ".", model.ext)
## A quick fix because names weren't used consistently.
model.names[model.names == "int.sep"] <- "int.psi"
## Number of models.
n.models <- length(model.ext)

## A data frame indicating which effects appear in which models.
model.df <- data.frame(
    fixed = c(rep(c("SST", "Periodic"), 8), rep(c("None", "SST", "Periodic"), 2)),
    spatially.varying = c(rep(rep(c("None", "SST"), each = 2), 4), rep(c("None", "SST"), each = 3)),
    spatiotemporal = c(rep(c("None", "Space", "Time"), each = 4),
                       rep("Space+Time", 4), rep("Space*Time", 6)))
              
## Creating a list, where each component is for a species.
fit <- rep.summary <- rand.summary <- d.full <- vector(mode = "list", length = n.species)
## Each component is itself a list with a component for each model.
for (i in 1:n.species){
    fit[[i]] <- rep.summary[[i]] <- rand.summary[[i]] <-
        d.full[[i]] <- vector(mode = "list", length = n.models)
    names(fit[[i]]) <- names(rep.summary[[i]]) <-
        names(rand.summary[[i]]) <- names(d.full[[i]]) <- paste("Model", 1:22)
}

## Extracting fixed-effect summaries, random effect summaries, and
## estimated spatiotemporal fields for each model/species
## combination. Some warnings are due to non-identifiable models that
## arise from attempting to fit complicated effects to species with
## low numbers of detections, and these are removed from consideration
## by the mdoel selection procedure.
for (i in 1:n.species){
    for (j in 1:n.models){
        filename <- paste("fits", "/fit-", model.ext[j], "-species-", i, ".RData", sep = "")
        if (file.exists(filename)){
            load(filename)
            fit[[i]][[j]] <- get(paste("fit.", model.names[j], sep = ""))
            rep.summary[[i]][[j]] <- summary(get(paste("sdrep.", model.names[j], sep = "")),
                                             type = "report")
            rand.summary[[i]][[j]] <- summary(get(paste("sdrep.", model.names[j], sep = "")),
                                              type = "random")
            d.full[[i]][[j]] <- get(paste("d.full.", model.names[j], sep = ""))
        }
    }
}

save(fit, rep.summary, rand.summary, d.full, model.df, file = "all-species-output.RData")
