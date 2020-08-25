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


fit <- vector(mode = "list", length = n.models)
rep.summary <- vector(mode = "list", length = n.models)
rand.summary <- vector(mode = "list", length = n.models)
d.full <- vector(mode = "list", length = n.models)
names(fit) <- names(rep.summary) <- names(rand.summary) <- names(d.full) <- model.ext

for (i in 1:n.models){
    cat(i)
    filename <- paste("fit-", model.names[i], ".RData", sep = "")
    load(filename)
    fit[[i]] <- get(paste("fit.", model.ext[i], sep = ""))
    rep.summary[[i]] <- summary(get(paste("sdrep.", model.ext[i], sep = "")), type = "report")
    rand.summary[[i]] <- summary(get(paste("sdrep.", model.ext[i], sep = "")), type = "random")
    d.full[[i]] <- get(paste("d.full.", model.ext[i], sep = ""))
}

out.name <- paste0("all-output", "-new"[new], ".RData")
save(fit, rep.summary, rand.summary, d.full, file = out.name)
