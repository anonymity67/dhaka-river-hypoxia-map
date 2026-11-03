# Dhaka river hypoxia risk map

Interactive map of the modelled probability that dissolved oxygen (DO) falls below 2 mg/L along Dhaka's peripheral rivers, by season.

- Data: 1,495 water-quality samples at 60 stations on 12 rivers, 2017-2021 (IEEE ClimateChain / Climate++ dataset).
- Model: Bayesian Bernoulli model fitted by INLA-SPDE (season effect, year effect, Matern spatial field). Leave-one-station-out AUC 0.91.
- Predictions are for a typical year, on OpenStreetMap river, canal, stream and drain lines within 5 km of a station.

## Files

- `index.html` - the page (Leaflet, no build step)
- `data/rivers.geojson` - river lines with mean predicted probability per season (`w` winter, `p` pre-monsoon, `m` monsoon, `o` post-monsoon), river group `g`, waterway type `t`, OSM name `n`
- `data/stations.json` - station name, river group, location, share of samples that were hypoxic (`h`)

## Limits

Hypoxia is a habitat condition, not a measurement of aquatic life. Sampling was about once per season. Values between stations are interpolated. This is research output, not a regulatory assessment.

## Data source and citation

Water-quality samples and station locations come from:

Sinha, P., Debnath, N., & Bhattacharjee, P. (2026). *Water Quality and Hydroclimatic Dataset for Rivers and Streams of Dhaka, Bangladesh* (Version 1) [Data set]. Mendeley Data. https://doi.org/10.17632/xftpv7jvx3.1 (CC BY 4.0).

The values in `data/` are derived from it by a statistical model (not the original measurements).

## Licences

- Code (`index.html`): MIT, see `LICENSE`.
- Derived data in `data/`: CC BY 4.0, see `DATA_LICENSE.md`.
- River lines and basemap: (c) OpenStreetMap contributors, ODbL.

## How to cite this map

Antar, A. (2026). Dhaka river hypoxia risk map [Interactive map]. GitHub. (Add the repository URL once published; add the paper reference when available.)

## Publish on GitHub Pages

1. Create a new public repository and upload these files.
2. Settings > Pages > Deploy from a branch > `main` / root.
3. The page appears at `https://<username>.github.io/<repository>/`.
