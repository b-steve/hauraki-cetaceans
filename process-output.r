model.names <- c("fixed", "fixed-p", "st", "st-p", "int", "int-psi", "int-p")
n.models <- length(model.names)

model.ext <- chartr("-", ".", model.names)

fit <- vector(mode = "list", length = n.models)
rep.summary <- vector(mode = "list", length = n.models)
rand.summary <- vector(mode = "list", length = n.models)
d.full <- vector(mode = "list", length = n.models)
names(fit) <- names(rep.summary) <- names(rand.summary) <- names(d.full) <- model.ext

for (i in 1:n.models){
    filename <- paste("fit-", model.names[i], ".RData", sep = "")
    load(filename)
    fit[[i]] <- get(paste("fit.", model.ext[i], sep = ""))
    rep.summary[[i]] <- summary(get(paste("sdrep.", model.ext[i], sep = "")), type = "report")
    rand.summary[[i]] <- summary(get(paste("sdrep.", model.ext[i], sep = "")), type = "random")
    d.full[[i]] <- get(paste("d.full.", model.ext[i], sep = ""))
}

save(fit, rep.summary, rand.summary, d.full, file = "all-output.RData")
