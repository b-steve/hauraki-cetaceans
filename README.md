# Redistribution of cetaceans based on ocean temperature using spatiotemporal mixed effects models

This repository contains data and code to accompany the following manuscript:

* Colbert, J., Stevenson, B. C., Bowen, M., and Constantine, R. C. (2026) Redistribution of cetaceans based on ocean temperature using spatiotemporal regression models. *Diversity and Distributions*, *32*(2), e70145. [(link)](https://doi.org/10.1111/ddi.70145)

## Data

Sighting data (`sighting.df`) and monthly temperatures (`month.temp`) are available in the file `sighting.RData`.

We split the Hauraki Gulf into 1381 discrete cells. The data frame `sighting.df` includes a row for each cell for every month from August 2000 to June 2019. The columns are as follows:
* `year` and `month`: The year and month corresponding to the row.
* `pixel_id`: An ID number for the cell corresponding to the row.
* `total.trips`: The total number of trips conducted in that month.
* Columns beginning with `n.trips`: The total number of trips that sighted a particular species in that month and cell. For example, `n.trips.cdolp` is the number of trips that sighted common dolphin.
* `av.vesselprob`: The cell's visitation probability, calculated by taking a random sample of trips and computing the proportion that visited the cell corresponding to the row.

The spatial coordinates of each of the 1381 cells are available in the
file `pixelcoord.RData`.

## Dependencies

The following R packages are required to run the code in this repository: `fields`, `INLA`, `RColorBrewer`, `TMB`, `sf`, and `sp`.

## Code

Code to fit the models described in the manuscript is available in `fit.r`. The models involve Gaussian random fields containing spatiotemporal random effects. Computation of model likelihoods requires approximating high-dimensional integrals, achieved using the R package `TMB`. The `TMB` model-fitting code is written in C++, and is available in `binomial_fit.cpp`, which is compiled and then executed for each fitted model within `fit.r`.

Such models can be computationally expensive to fit, and ours require more RAM than is typically available on a standard desktop or laptop computer. Our code is structured for use on the high-performance computing resources provided by New Zealand's eScience Infrastructure. Fitting all models we consider to a species' sightings can be achieved by running the following at a command line:
```bash
R --vanilla < fit.r --args i
```
Replacing `i` with the following values to select the species:
- 1: Bryde's whale
- 2: Common dolphin
- 3: Bottlenose dolphin
- 4: Killer whale
- 5: Whale
- 6: Bryde's + whale

To save the bother of anyone else having to fit these models, the resulting `.RData` files (one for each combination of model and species) is available in the `fits/` directory. The code in `process-model-files.r` processes the model files in the `fit/` directory, combining all the necessary information into a single file, `all-species-output.RData`. This file is not available in this repositry because it is quite large (approx 300 MB), so it must be created locally by the user.

## Inference

Code to extract the inference available from our models is contained
in `model-output.r`. Model output includes
- model selection,
- plotting estimates of spatial effects, including spatially varying effects of temperature on cetacean occurrence,
- plotting estimates of temporal effects,
- animated GIFs of estimated spatiotemporal estimates of species occurrence,
- plotting spatial variation in effort, and
- plotting estimated spatial and temporal covariance functions.


