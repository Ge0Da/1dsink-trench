# Requirements

## MATLAB Version
- MATLAB R2022a or later

## Required Toolboxes
- Mapping Toolbox (GeoTIFF reading, coordinate transforms)
- Signal Processing Toolbox (spectral analysis)
- Statistics and Machine Learning Toolbox (Spearman correlation, bootstrap)

## Input Data
- `data/kamchatka_gmrt.tif` — GMRT bathymetric DEM (~4.7 GB, not included in repo)
  Download from: https://www.gmrt.org/ (Kuril-Kamchatka region, ~41 m resolution)
- `data/kamchatka_samples.xlsx` — geochemical sample table (included)

## Verification
Run `reproduce_all.m` — it checks toolbox availability before starting.
