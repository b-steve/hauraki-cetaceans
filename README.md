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
