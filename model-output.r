## Loading in packages.
library(TMB)
library(INLA)
library(RColorBrewer)
library(fields)
library(sp)
library(sf)

## Loading in data and helper functions.
source("helper-funs.r")
## Setting up sighting counts for the different species categories.
y <- as.matrix(sighting.df[, c(8:11, 15)])
## The final column is the Bryde's plus general whale categories.
y <- cbind(y, y[, 1] + y[, 5])

## The following data frame indicates which effects appear in which
## models.
model.df

## AICs for model selection.
aics <- calc.aics()
## This is a list with a bunch of components:
## - $aic is the AIC value for each species-model pair.
## - $diff is the difference between a model's AIC and the model with
##   the best AIC for that species, and reflects the AIC table in the
##   manuscript.
## - $best has a TRUE entry in the best model (by AIC) for each
##   species.
aics

## For lots of the codes below, use these values for the species argument:
## 1 = "byrde",
## 2 = "cdolp",
## 3 = "bdolp",
## 4 = "orca",
## 5 = "whale"
## 6 = "brydeplus"

## And choose a model number from model.df

## Plotting estimated relative density for a particular month/species/model.
plot.surf(species = 1, model = 21, month = 123)
## Plotting estimated spatially varying effect of temperature for a
## particular species/model. Only makes sense for models that include
## such an effect, where the 'spatially.varying' column is "SST" in
## model.df.
plot.surf(species = 1, model = 16, surf = "int")
## Plotting estimated spatial effect. This only makes sense for models
## that have a single constant wiggly spatial field over time, where
## the 'spatiotemporal' column is "Space" or "Space+Time" in model.df.
plot.surf(species = 1, model = 16, surf = "omega")

## For example, here's a plot of the spatially varying effect of
## temperature for all species.
par(mfrow = c(3, 2))
## Just using the model 16 for these, which isn't actually the best
## model for every species.
for (i in 1:6){
   plot.surf(species = i, model = 16, surf = "int")
}

## Here's a plot we used for the paper for the species that exhibited
## spatially varying effects.
models <- apply(aics$best, 1, which)
par.save <- par(mfrow = c(2, 2), mar = c(3, 4, 2, 2), oma = c(0.1, 0.1, 0.1, 1.1))
plot.title <- c("A. Bryde's whale", "B. Common dolphin", "C. Whale", "D. Bryde's whale +")
k <- 1
for (i in c(1, 2, 5, 6)){
    plot.surf(species = i, model = models[i], surf = "int", main = plot.title[k])
    k <- k + 1
}
par(par.save)

## Plotting estimated temporal effect. This only makes sense for
## models that have a wiggly temporal process over time, where the
## 'spatiotemporal' column is "Time" or "Space+Time" in model.df.
plot.temporal(species = 1, model = 16)

## Make a distribution gif for a species/model combination. The code
## below creates a gif for each species, using the best model for that
## species. Takes a while to run. Note that this function requires the
## command-line tool convert. You can check by running the following
## in R:
## system("convert --version")
## If available, you will see output related to your version of
## convert. The gif is saved to the current R working directory.
for (i in 1:6) save.gif(species = i, model = which(aics$best[i, ]))

## Plot survey effort.
plot.effort()

## Plotting changes in occurrence probability at a specific locations
## over time. First, select location IDs from the plot below to
## include in the plot of temporal trends.
loc.id.selector()
## The locations IDs we've selecte from above are here.
loc.ids <- c(310, 311, 292, 295)
## A vector indicating which model we plot for each species. Here we
## use whichever is best by AIC.
models <- apply(aics$best, 1, which)
## Plots of temporal trends at the selected locations.
plot.temporal.at.locs(loc.ids, models)

## Plotting spatial covariance functions. Selected models must include
## spatial random effects.
plot.spatial.cov(models)

## Plotting temporal covariance functions. Selected models must inclue
## temporal random effects.
plot.temporal.cov(models)


