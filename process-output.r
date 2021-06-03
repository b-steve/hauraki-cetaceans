new <- TRUE
if (new){
    model.names <- c("cf", "fixed", "int", "st")
} else {
    model.names <- c("fixed", "fixed-p", "st", "st-p", "int", "int-psi", "int-p")
}
n.models <- length(model.names)

model.ext <- chartr("-", ".", model.names)
if (new){
    model.names <- paste0(model.names, "-new")
}


library(TMB)
smalltri <- TRUE
## Number of species.
n.species <- 6
## Filename indicator for each of the models.
model.ext <- c("cf", "fixed-p", "fixed", "int-p", "int-sep", "int", "st-p", "st")
## R object indicator for each of the models.
model.names <- chartr("-", ".", model.ext)
## Because I didn't name things consitently...
model.names[5] <- "int.psi"
## Number of models.
n.models <- length(model.ext)

## Creating a list, where each component is for a species.
fit <- rep.summary <- rand.summary <- d.full <- vector(mode = "list", length = n.species)
## Each component is itself a list with a component for each model.
for (i in 1:n.species){
    fit[[i]] <- rep.summary[[i]] <- rand.summary[[i]] <-
        d.full[[i]] <- vector(mode = "list", length = n.models)
    names(fit[[i]]) <- names(rep.summary[[i]]) <-
        names(rand.summary[[i]]) <- names(d.full[[i]]) <- model.ext
}

for (i in 1:n.species){
    for (j in 1:n.models){
        cat(i, j, "\n")
        current.obj <- ls()
        filename <- paste("fits-newest/fit-", model.ext[j], "-species-", i, ".RData", sep = "")
        if (file.exists(filename)){
            load(filename)
            fit[[i]][[j]] <- get(paste("fit.", model.names[j], sep = ""))
            rep.summary[[i]][[j]] <- summary(get(paste("sdrep.", model.names[j], sep = "")),
                                             type = "report")
            rand.summary[[i]][[j]] <- summary(get(paste("sdrep.", model.names[j], sep = "")),
                                              type = "random")
            d.full[[i]][[j]] <- get(paste("d.full.", model.names[j], sep = ""))
        }
        new.obj <- ls()
        rm(list = new.obj[!new.obj %in% current.obj])
    }
}

out.name <- paste0("all-species-output", "-smalltri"[smalltri], ".RData")
save(fit, rep.summary, rand.summary, d.full, file = out.name)
