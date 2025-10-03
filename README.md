# Redistribution of cetaceans based on ocean temperature using spatiotemporal mixed effects models

This repository contains data and code to accompany the following manuscript:

* Colbert, J., Stevenson, B. C., Bowen, M., and Constantine, R. C. (in submission) Redistribution of cetaceans based on ocean temperature using spatiotemporal mixed effects models.

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

## Code

Code to fit the models described in the manuscript is available in `fit.r`. The models involve Gaussian random fields containing spatiotemporal random effects. Computation of model likelihoods requires approximating high-dimensional integrals, achieved using the R package `TMB`. 

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

