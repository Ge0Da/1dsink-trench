# 1DSink Along-Trench Morphometric Framework

One-dimensional along-trench morphometric decomposition and sediment-trapping analysis of the Kuril-Kamchatka Trench (~3,100 km), validated against independent surface-sediment geochemistry. The pipeline extracts the trench axis from bathymetry, decomposes the along-axis depth profile into wavelength bands, segments it into morphological domains, and combines low slope, roughness, and concavity into a composite Trapping Index (TI).

## Citation

If you use this code, please cite the associated publication and the archived software release.

Fernández-Blanco, D., Dasgupta, S., et al. Morphology segregates contaminant accumulation along hadal trenches (under review).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Quick start

```matlab
cd <repository root>
reproduce_all          % full pipeline + supplementary analyses; checks toolboxes first
```

Or run just the main figure pipeline:

```matlab
addpath('utils')
run_master_figure
```

Outputs (written to `output/`): the 5-panel master figure, a QC diagnostic figure, the band-decomposition CSV, the projected sample positions, and a `.mat` file holding all pipeline results.

## Structure

```
.
├── config_master_figure.m     Configuration: paths, spectral bands, thresholds
├── run_master_figure.m        Driver script (full pipeline)
├── reproduce_all.m            Reproduce all results from scratch (checks toolboxes)
├── LICENSE.txt                MIT License
├── REQUIREMENTS.md            MATLAB version, toolboxes, input data
│
├── analysis/                  Supplementary analyses (run after the main pipeline)
│   ├── calibrate_TI_weights.m TI weight calibration (simplex grid search + LOOCV)
│   ├── geochem_crosscheck.m   Spearman correlations: TI vs contaminant concentrations
│   ├── parameter_sweep.m      Constrained segmentation sensitivity analysis
│   ├── stationarity_test.m    Scale-invariance check of domain boundaries
│   └── test_convergence.m     Convergence-scoring verification
│
├── utils/                     Functions (24, each documented in its header)
│   ├── autocorrDomainCount.m  Data-driven minimum domain count from the ACF
│   ├── calib_config.m         Calibration-target metal list
│   ├── computeArcLength.m     Along-axis arc length
│   ├── computeFlowRouting.m   Flow direction / accumulation on the inverted DEM
│   ├── computeSinkMetrics.m   Sliding-window roughness + Trapping Index
│   ├── computeTrenchBands.m   Multi-band Gaussian decomposition
│   ├── convergenceScore.m     Multi-proxy boundary convergence (contrast ratios)
│   ├── drawSegmentLines.m     Domain-boundary line rendering
│   ├── drawTIcurve.m          Composite Trapping-Index banner
│   ├── extractMainDivide.m    Drainage-divide extraction
│   ├── extractTrenchAxis.m    Trench-axis extraction pipeline
│   ├── orderPointsByConnectivity.m  PCA-based ordering of axis points
│   ├── plotSampleLabels.m     Non-overlapping sample labels
│   ├── plotSinkMetricsDiagnostic.m  QC diagnostic figure (6 panels)
│   ├── plotTrenchBands_1DSink.m     Master figure (5 panels)
│   ├── projectSamplesOnAxis.m Project samples onto the axis (along/lateral)
│   ├── refineToLocalMinimum.m Snap the axis to full-resolution depth minima
│   ├── reprojectSamplesToGrid.m  Reproject sample coordinates to the grid CRS
│   ├── resampleToCoarse.m     DEM resampling to coarse resolution
│   ├── sampleFromDEM.m        Depth sampling from the DEM
│   ├── saveFigure.m           Multi-format figure export (PNG + SVG)
│   ├── segmentDomains.m       Curvature-based domain segmentation
│   ├── targetMean.m           Scaled mean of the calibration-target metals
│   └── trimBands.m            Display-range trimming of bands
│
├── data/                      Input data
│   └── kamchatka_samples.xlsx Geochemical sample table (17 stations)
│       (sediment_texture_TOC.csv — per-station mud% and TOC)
│       (kamchatka_gmrt.tif — GMRT bathymetric DEM, ~4.7 GB, NOT included; see below)
│
└── output/                    Generated products (figures, band CSV, results .mat)
```

The bathymetric DEM (`data/kamchatka_gmrt.tif`, ~4.7 GB) is not distributed with this repository; download it from GMRT (see Data).

## Pipeline

1. **Trench-axis extraction** from inverted bathymetry via a dual-resolution drainage-divide method (500 m coarse, then refined to local depth minima on the native DEM).
2. **Sample projection** onto the trench axis (17 geochemical stations).
3. **Multi-band Gaussian decomposition** separating the structural envelope (>10 km), basin-scale relief (2-20 km), and residual roughness (<2 km).
4. **Curvature-based domain segmentation** with a stability audit: sign changes in the D_{2-20} curvature band partition the axis into concave (trap) and convex (bypass) domains. Boundary stability is tested under +/-30% parameter perturbation.
5. **Multi-proxy convergence scoring**: each domain boundary is scored across curvature, slope gradient, and roughness via contrast ratios in +/-20 km windows. The composite strength (geometric mean of the three ratios) measures multi-scale morphological support.
6. **Composite Trapping Index**: TI combines normalised low slope, roughness, and concavity. The operational weights (0.47, 0.40, 0.13) are calibrated against the four-element heavy-metal mean (Hg, As, Cd, Pb; n = 12); see `analysis/calibrate_TI_weights.m` and `analysis/geochem_crosscheck.m`.
7. **Figure generation** with boundary line weight and opacity scaled by composite convergence strength.

## Supplementary analyses

Scripts in `analysis/` resolve the project root automatically:

```matlab
addpath('utils')
run('analysis/calibrate_TI_weights')   % TI weight optimisation + LOOCV
run('analysis/geochem_crosscheck')     % Spearman correlations vs contaminants
run('analysis/parameter_sweep')        % Constrained segmentation sensitivity
run('analysis/stationarity_test')      % Domain-boundary scale-invariance
run('analysis/test_convergence')       % Convergence-score verification
```

## Requirements

MATLAB R2022a or later. Toolboxes: Mapping, Signal Processing, Statistics and Machine Learning. See `REQUIREMENTS.md`. `reproduce_all.m` checks toolbox availability before starting.

## Data

- **Bathymetry**: GMRT synthesis (https://www.gmrt.org/), Kuril-Kamchatka region, gridded at ~41 m. Export the region as GeoTIFF to `data/kamchatka_gmrt.tif`.
- **Geochemistry**: `data/kamchatka_samples.xlsx` (17 stations) and `data/sediment_texture_TOC.csv` (per-station mud% and TOC), both included.

## Key parameters

Tuneable parameters live in `config_master_figure.m`, each annotated with a sensitivity rating.

| Parameter | Value | Sensitivity | Notes |
|---|---|---|---|
| Spectral bands | 2, 5, 10, 20, 50, 200, 1000 km | high | Controls all downstream decomposition |
| Window size | 50 km | moderate | Roughness computation window |
| Edge exclusion | 150 km | moderate | Quarantines filter ramp-up artefacts |
| TI weights | 0.47 : 0.40 : 0.13 | moderate | Calibrated against heavy-metal data |
| Deadband multiplier | 1.5 x MAD | moderate | Curvature ambiguity threshold |
| Convergence window | +/-20 km | low | Matches the D_{2-20} maximum wavelength |

## License

MIT License (see `LICENSE.txt`). Copyright (c) 2026 David Fernández-Blanco, IDSSE-CAS.

## Contact

David Fernández-Blanco — davidfb@idsse.ac.cn
Deep Sea Geophysics and Resources Research Laboratory, IDSSE-CAS, Sanya, China
