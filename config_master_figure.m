function cfg = config_master_figure()
%CONFIG_MASTER_FIGURE  Central configuration for 1DSink Master Figure production
%
%   cfg = config_master_figure()
%
%   Returns a configuration structure containing tunable parameters
%   for the 1DSink Along Strike Master Figure and QC Diagnostic pipelines.
%   Edit values below to adapt the analysis to the trench system / dataset.
%
%   ORGANISATION:
%     Section 1 — File paths (DEM, samples, output)
%     Section 2 — Axis extraction (coarse resolution, smoothing, trimming)
%     Section 3 — Signal decomposition (resampling, filters, wavelengths)
%     Section 4 — Sink metrics & Trapping Index (windowing, weights, β)
%     Section 5 — Domain segmentation (deadband, persistence, confidence)
%     Section 6 — Sample handling (column names, projection)
%     Section 7 — Figure formatting (sizes, fonts, colours, line widths)
%     Section 8 — Output (formats, resolution, verbosity)
%     Appendix  — Internal constants (documented, not tunable)
%
%   DIFFICULTY INDICATORS:
%     ★ Routine — safe to change for any dataset
%     ★★ Intermediate — understand the downstream effect first
%     ★★★ Advanced — changes affect the entire spectral architecture
%
%   USAGE:
%     1. Edit the values in this file.
%     2. Run:  run_master_figure
%     The driver reads this config and passes values to every function.
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date:   2026-03-27
%   Version: 1.1

%% ========================================================================
%  SECTION 1: FILE PATHS  ★
%  ========================================================================
%  Paths to input data and output directory. Relative paths are resolved
%  from the directory containing run_master_figure.m.

% -------------------------------------------------------------------------
% DEM file (GeoTIFF or NetCDF readable by TopoToolbox GRIDobj)
% This must cover the trench and its surrounding slopes.
% Typical sources: GEBCO_2025, SRTM15+, regional multibeam compilations.
% -------------------------------------------------------------------------
cfg.demFile = fullfile('data', 'kamchatka_gmrt.tif');

% -------------------------------------------------------------------------
% Sample file (CSV or Excel table)
%   REQUIRED columns: coordinates (see Section 6 for column names)
%   REQUIRED: at least one name/ID column
%   Samples are mandatory — errors if this file is missing.
% -------------------------------------------------------------------------
cfg.sampleFile = fullfile('data', 'kamchatka_samples.xlsx');

% -------------------------------------------------------------------------
% Output directory — created automatically if it does not exist.
% All figures, MAT files, and CSVs are saved here.
% -------------------------------------------------------------------------
cfg.outputDir = 'output';


%% ========================================================================
%  SECTION 2: AXIS EXTRACTION
%  ========================================================================
%  Parameters controlling the drainage-divide method for trench axis
%  extraction (dual-resolution approach: coarse extraction + fine
%  refinement).

% -------------------------------------------------------------------------
% ★  Coarse-phase DEM resolution (metres)
%
%  WHAT IT DOES:
%    The trench axis is first extracted on a coarsened DEM for speed, then
%    refined to the full-resolution DEM. This sets the coarsened cell size.
%
%  PHYSICAL MEANING:
%    At 500 m, features smaller than ~1 km are suppressed, allowing the
%    algorithm to capture the trench-scale (100+ km) drainage divide
%    without being confused by small-scale roughness (horst–graben fabric,
%    seamounts).
%
%  IF TOO SMALL (< 300 m):
%    The divide extraction sees too much topographic complexity, potentially
%    fragmenting the divide or routing it along local ridges. Computation
%    time increases quadratically (250 m grid has 4x more cells than 500 m).
%
%  IF TOO LARGE (> 1500 m):
%    Narrow trench features and forearc ridges are smoothed away. The
%    algorithm may fail to identify the correct divide, especially in
%    very narrow trenches (<~5 km).
%
%  COMMON RANGE: 200–1000 m
%  DEFAULT: 500 m
% -------------------------------------------------------------------------
cfg.axis.coarseRes_m = 500;

% -------------------------------------------------------------------------
% ★  Minimum drainage area for stream network (km²)
%
%  WHAT IT DOES:
%    Sets the flow-accumulation threshold for the coarse-resolution stream
%    network. Only channels draining ≥ this area are used for divide
%    extraction.
%
%  IF TOO SMALL (< 30 km²):
%    Dense network; divide algorithm may find internal divides within the
%    trench-ward drainage rather than the main bounding divide.
%
%  IF TOO LARGE (> 500 km²):
%    Too few channels to define a meaningful divide. Algorithm may fail.
%
%  INTERACTION: Effective cell count = minAreaKm2 × 1e6 / coarseRes_m².
%    At 500 m, 100 km² = 400 cells. Keep effective count ≥ 50 cells.
%
%  COMMON RANGE: 10–1000 km²
%  DEFAULT: 100 km²
% -------------------------------------------------------------------------
cfg.axis.minAreaKm2 = 100;

% -------------------------------------------------------------------------
% ★  Local refinement radius (pixels at full DEM resolution)
%
%  WHAT IT DOES:
%    After coarse-phase extraction, each axis point is snapped to the
%    nearest local depth minimum within this radius on the full-res DEM.
%
%  IF TOO SMALL (1):
%    May fail to find the true thalweg if coarse estimate is offset.
%
%  IF TOO LARGE (> 7):
%    May "jump" to adjacent features (nearby channel, slope basin).
%
%  COMMON RANGE: 1–10 pixels
%  DEFAULT: 3  (±1.5 km at 500 m resolution)
% -------------------------------------------------------------------------
cfg.axis.localRefineRadius_cells = 3;

% -------------------------------------------------------------------------
% ★  Axis planform smoothing window (number of points)
%
%  WHAT IT DOES:
%    Removes pixel-scale zigzag artefacts from the grid-based extraction.
%    Applied as a moving average to (x, y) coordinates, then depth is
%    re-sampled from the DEM at smoothed locations.
%
%  EFFECTIVE WAVELENGTH: smoothWindow × point_spacing.
%    At 100 m spacing, 11 points = 1.1 km.
%
%  IF TOO SMALL (3): Staircase artefacts remain, creating spurious
%    high-frequency noise in derivatives.
%  IF TOO LARGE (> 31): Real along-trench bends are smoothed out.
%
%  NOTE: Must be ODD (even windows are asymmetric). The driver enforces
%  this by rounding up to the next odd integer if necessary.
%
%  COMMON RANGE: 3–51
%  DEFAULT: 11
% -------------------------------------------------------------------------
cfg.axis.smoothWindow_samples = 11;

% -------------------------------------------------------------------------
% ★  Axis trimming (km along axis)
%
%  WHAT IT DOES:
%    Crops the extracted axis to a specified along-axis range. For final
%    output after exploratory work, record the known-good range here.
%
%  FORMAT:
%    []             — Use the full extracted axis (no trimming)
%    [start, end]   — Keep only the range from start_km to end_km
%
%  EXAMPLE: cfg.axis.axisTrim_km = [200, 1800];
%    Keeps only km 200 to km 1800. Arc length is recomputed from zero.
%
%  HOW TO DETERMINE: Run the full-axis pipeline at least once. Record the
%  area of interest selected range in km. Enter it here.
%
%  DEFAULT: [] (use full axis)
% -------------------------------------------------------------------------
cfg.axis.axisTrim_km = [];

% -------------------------------------------------------------------------
% ★  Refine to full DEM resolution (logical)
%
%  If true, the coarse-phase axis is refined by snapping each point to the
%  local depth minimum on the full-resolution DEM. Set false only if the
%  DEM has a single resolution (no coarsening step was applied).
%
%  DEFAULT: true
% -------------------------------------------------------------------------
cfg.axis.refineToFullRes = true;


%% ========================================================================
%  SECTION 3: SIGNAL DECOMPOSITION
%  ========================================================================
%  Parameters controlling the multi-band Gaussian decomposition of the
%  along-trench depth profile. Changes here affect ALL downstream
%  product: band-pass amplitudes, derivatives, metrics, domains, and
%  Trapping Index.

% -------------------------------------------------------------------------
% ★★  Uniform resampling interval (metres)
%
%  WHAT IT DOES:
%    The extracted axis has irregular point spacing. Before spectral
%    decomposition, the depth profile is resampled onto a uniform grid.
%
%  NYQUIST CONSTRAINT: minimum resolvable wavelength = 2 × spacing.
%    At 100 m, features ≥ 200 m are resolvable.
%
%  RULE: Match to DEM resolution. Never go below the DEM native resolution
%    (creates phantom detail through interpolation).
%
%  COMMON RANGE: 50–1000 m
%  DEFAULT: 100 m  (for ~100 m multibeam DEMs)
% -------------------------------------------------------------------------
cfg.signal.resampleSpacing_m = 100;

% -------------------------------------------------------------------------
% ★  Interpolation method for resampling
%
%  'makima'  — Modified Akima: shape-preserving, no overshoot (DEFAULT)
%  'pchip'   — Piecewise cubic Hermite: smooth, minimal overshoot
%  'linear'  — Fast but introduces corners → spurious high-frequency content
%
%  Use 'makima' unless MATLAB is older than R2019a (then 'pchip').
%  Avoid 'linear' for spectral analysis.
%
%  DEFAULT: 'makima'
% -------------------------------------------------------------------------
cfg.signal.interpMethod = 'makima';

% -------------------------------------------------------------------------
% ★★  Edge artefact exclusion zone (km)
%
%  WHAT IT DOES:
%    Gaussian filters produce artefacts at signal edges. This parameter
%    defines a zone at each end of the profile flagged as unreliable via
%    an edgeMask (false at edges, true in interior).
%
%  RULE OF THUMB: edgeTrim >= 2.5 x sigma_max, where sigma_max = lambda_max / (2*pi)
%    from the longest wavelength used in the SEGMENTATION bands (not display bands).
%    For lambda = 200 km (L200): sigma = 31.8 km, half-width = 127 km -> 150 km is safe.
%    For lambda = 1000 km (L1000): sigma = 159 km, half-width = 637 km -> 650 km would
%      be needed, but this is impractical for most profiles. Accept that the
%      mega-trend band (LP>1000) has edge leakage; it affects Panel A display
%      only, not segmentation or TI.
%
%  IF TOO SMALL (< 50 km): Edge artefacts in D2-20 segmentation band leak
%    into domain classification near profile ends.
%  IF TOO LARGE (> 300 km): Excessive data loss (problematic for trenches
%    shorter than ~1000 km).
%
%  COMMON RANGE: 50-200 km
%  DEFAULT: 150 km (safe for segmentation and metric bands up to L200)
% -------------------------------------------------------------------------
cfg.signal.edgeTrim_km = 150;

% -------------------------------------------------------------------------
% ★★  Derivative smoothing window (km)
%
%  WHAT IT DOES:
%    Numerical derivatives amplify high-frequency noise. This sets the
%    width of a moving-average smoother applied to first and second
%    derivatives of band-pass components.
%
%  EFFECTIVE WINDOW: derivWindow_km / spacing.
%    At 100 m and 2 km → 20-point moving average.
%
%  IF TOO SMALL (< 0.5 km): Noisy derivatives, chaotic curvature.
%  IF TOO LARGE (> 5 km): Blurs real 5–10 km features.
%
%  NOTE: Segmentation curvature uses an additional 3× multiplier
%    internally (hardcoded — see Appendix).
%
%  COMMON RANGE: 0.5–10 km
%  DEFAULT: 2 km
% -------------------------------------------------------------------------
cfg.signal.derivWindow_km = 3;

% -------------------------------------------------------------------------
% ★★★  Filter wavelength cutoffs (km) — THE SPECTRAL ARCHITECTURE
%
%  WHAT IT DOES:
%    Defines the Gaussian low-pass filter cutoff wavelengths that form the
%    multi-band decomposition. Band-pass signals are formed by subtraction.
%
%  BAND STRUCTURE (from these cutoffs):
%    D_{<10}  = raw − LP_10      S2S-scale roughness
%    D_{2-10} = LP_2 − LP_10     PRIMARY slope band
%    D_{5-10} = LP_5 − LP_10     PRIMARY curvature band
%    D_{2-20} = LP_2 − LP_20     Segmentation band (domain boundaries)
%    D_{10-50} = LP_10 − LP_50   Tectonic band 1 (flexural)
%    D_{50-200} = LP_50 − LP_200 Tectonic band 2 (regional)
%    D_{>1000} = LP_1000         Mega-scale trend
%
%  PHYSICAL BASIS FOR DEFAULT BOUNDARIES:
%    10 km — Typical trench width (5–15 km; Stern, 2002). Features < 10 km
%            are sedimentary architecture. Features > 10 km are plate-scale.
%    50 km — Typical spacing of forearc structural segments and outer-rise
%            fault blocks (Kobayashi et al., 1998).
%   200 km — Characteristic wavelength of along-strike slab dip variations
%            (Schellart et al., 2007).
%
%  WHEN TO MODIFY:
%    Narrow trenches (< 5 km width): reduce L10 to 5 km.
%    Wide trenches with thick fill (> 20 km): increase L10 to 15–20 km.
%
%  ⚠ CAUTION: Changing ANY wavelength affects ALL downstream products.
%    Re-validate the entire analysis after any change.
%
%  DEFAULT: struct('L1',1,'L2',2,'L5',5,'L10',10,'L20',20,'L50',50,
%                  'L200',200,'L1000',1000)
% -------------------------------------------------------------------------
cfg.signal.lambda_km = struct( ...
    'L1',    1, ...     % 1 km cutoff
    'L2',    2, ...     % Slope band lower bound
    'L5',    5, ...     % Curvature band lower bound
    'L10',  10, ...     % *Selected, tectonic boundary wavelength for KKT*
    'L20',  25, ...     % Segmentation upper bound
    'L50',  100, ...    % Tectonic band 1
    'L200', 200, ...    % Tectonic band 2
    'L1000', 1000);     % Regional flexure


%% ========================================================================
%  SECTION 4: SINK METRICS & TRAPPING INDEX
%  ========================================================================
%  Parameters for windowed roughness metrics (RMS, MAD, spectral slope β)
%  and the composite Trapping Index.

% -------------------------------------------------------------------------
% ★★  Moving window size (km)
%
%  WHAT IT DOES:
%    Controls the spatial scale at which roughness metrics and TI are
%    computed. Each metric is evaluated in a sliding window of this width.
%
%  PHYSICAL MEANING: A 50 km window captures the roughness character of
%    one forearc structural segment. At 100 m spacing, the window contains
%    ~500 points and ~5–25 cycles of the 2–10 km roughness.
%
%  IF TOO SMALL (< 25 km): Too few cycles for stable statistics.
%  IF TOO LARGE (> 100 km): Local variations are averaged out.
%
%  INTERACTION: Should be ≥ 2× minDomainLength_km.
%
%  COMMON RANGE: 20–200 km
%  DEFAULT: 50 km
% -------------------------------------------------------------------------
cfg.metrics.windowSize_km = 50;

% -------------------------------------------------------------------------
% ★  Window overlap fraction (0–1)
%
%  Higher overlap → smoother curves (but spatially autocorrelated).
%  Lower overlap → independent windows (better for hypothesis testing).
%
%  COMMON RANGE: 0.0–0.95
%  DEFAULT: 0.8  (smooth curves for figure production)
% -------------------------------------------------------------------------
cfg.metrics.windowOverlap = 0.85;

% -------------------------------------------------------------------------
% ★★  Trapping Index component weights
%
%  TI = (w_slope × LowSlope + w_rough × Roughness + w_curv × Concavity)
%       / (w_slope + w_rough + w_curv)
%
%  COMPONENT INTERPRETATION:
%    LowSlope:  Low axial gradient favours ponding.    Weight: w_slope
%    Roughness: High roughness = more mini-basins.     Weight: w_rough
%    Concavity: Concave-up = basin geometry.           Weight: w_curv
%
%  DEFAULT RATIONALE (1.4 : 1.2 : 0.4):
%    Calibrated against heavy-metal geochemistry (n=12, KKT).
%    Conservative shift from original 1.0:1.2:0.8 toward bootstrap
%    medians [0.55, 0.35, 0.10].  LOOCV rho = 0.455 (n=12) shows
%    optimised weights overfit; this is a halfway compromise.
%
%  MODIFICATIONS:
%    Canyon-dominated systems: increase w_slope
%    Ponded-basin systems: increase w_curv
%    Rough outer-rise systems: decrease w_rough
%    Set any weight to 0 to remove that component entirely.
%
%  COMMON RANGE: 0.0–5.0 each
%  DEFAULT: w_slope=1.0, w_rough=1.2, w_curv=0.8
% -------------------------------------------------------------------------
cfg.metrics.w_slope = 1.4; %Before 1.0
cfg.metrics.w_rough = 1.2; %Before 1.2
cfg.metrics.w_curv  = 0.4; %Before 0.8

% -------------------------------------------------------------------------
% ★★  TI normalization mode
%
%  'absolute'   — SPATIALLY STATIONARY (default).  TI at each window is
%                  determined solely by absolute physical thresholds below.
%                  Trimming the profile does NOT change any TI value.
%  'percentile' — non-stationary alternative.  TI depends on global percentiles
%                  across all windows.  Trimming changes TI everywhere.
%
%  DEFAULT: 'absolute'
% -------------------------------------------------------------------------
cfg.metrics.normMode = 'absolute';

% -------------------------------------------------------------------------
% ★★  Absolute TI thresholds (used when normMode = 'absolute')
%
%  These define physically meaningful scales for trench-floor processes.
%  A slope of 15 m/km is steep for a trench floor regardless of what
%  other slopes exist in the profile.
%
%  SLOPE (m/km): slopeLo → score 1 (flat = trapping)
%                slopeHi → score 0 (steep = bypass)
%    Typical trench floors: 0.5–20 m/km
%    DEFAULT: [0.5, 15.0]
%
%  ROUGHNESS (m): roughLo → score 0 (smooth)
%                 roughHi → score 1 (rough = mini-basins)
%    RMS of D_{2-10} band.  Typical range: 1–100 m.
%    DEFAULT: [2.0, 80.0]
%
%  CONCAVITY (m/km²): curvLo → score 0 (convex / neutral)
%                     curvHi → score 1 (concave-up = ponding)
%    Positive = concave-up in our convention (sign-flipped d2D).
%    DEFAULT: [-0.01, 0.05]
%
%  ⚠ TUNING: Run once with verbose=true to see raw metric ranges.
%    Adjust thresholds if your trench has very different scales.
% -------------------------------------------------------------------------
cfg.metrics.slopeLo_mpkm  = 0.5;
cfg.metrics.slopeHi_mpkm  = 15.0;
cfg.metrics.roughLo_m     = 2.0;
cfg.metrics.roughHi_m     = 80.0;
cfg.metrics.curvLo_mpkm2  = -0.01;
cfg.metrics.curvHi_mpkm2  = 0.05;

% -------------------------------------------------------------------------
% ★  Percentile-mode settings (used ONLY when normMode = 'percentile')
%
%  DEFAULT: pctileLow=5, pctileHigh=95
% -------------------------------------------------------------------------
cfg.metrics.pctileLow  = 5;
cfg.metrics.pctileHigh = 95;

% -------------------------------------------------------------------------
% ★★  Spectral slope (β) detrending method
%
%  'lp50'   — Remove wavelengths > 50 km. Isolates 2–50 km roughness.
%             (DEFAULT for standard 1DSink analysis)
%  'lp1000' — Remove wavelengths > 1000 km. Preserves broader spectrum
%             but more sensitive to flexural wavelengths.
%
%  ⚠ β must NEVER be computed on band-passed data (e.g., D_{2-10}).
%    Band-pass filtering imposes its own spectral roll-off.
%    Reference: Goff & Jordan (1988), JGR.
%
%  DEFAULT: 'lp50'
% -------------------------------------------------------------------------
cfg.metrics.betaDetrend = 'lp50';

% -------------------------------------------------------------------------
% ★★  Spectral slope (β) fitting band (cycles/km)
%
%  Frequency range for the power-law fit PSD(f) ∝ f^(−β).
%
%  CONVERSION: λ = 1/f.
%    [0.02, 0.5] cyc/km → [2 km, 50 km] wavelength.
%
%  IF BAND TOO NARROW: Few spectral points; β is unstable.
%  IF BAND TOO WIDE: Detrending artefacts contaminate the fit.
%
%  COMMON RANGE: [0.01–0.1, 0.3–1.0] cyc/km
%  DEFAULT: [0.02, 0.5]
% -------------------------------------------------------------------------
cfg.metrics.betaFitBand_cyc_km = [0.02, 0.5];


%% ========================================================================
%  SECTION 5: DOMAIN SEGMENTATION
%  ========================================================================
%  Parameters for the hysteresis-based curvature domain segmentation.

% -------------------------------------------------------------------------
% ★★  Curvature noise floor method (deadband mode)
%
%  Determines how the "deadband" — the zone around zero curvature where
%  sign is ambiguous — is computed.
%
%  'mad'        — Robust: uses Median Absolute Deviation of curvature,
%                 scaled by 1.4826 × deadbandSigmas. Adapts to noise
%                 level of the data. (DEFAULT, recommended)
%  'percentile' — Uses deadbandFactor × P80(|curvature|). Simpler but
%                 less robust.
%  'fixed'      — Uses deadbandAmplitude directly. Exact reproducibility
%                 independent of data statistics, but requires manual
%                 tuning per dataset.
%
%  DEFAULT: 'mad'
% -------------------------------------------------------------------------
cfg.segmentation.deadbandMode = 'mad';

% -------------------------------------------------------------------------
% ★  Fixed deadband amplitude (m/km²) — used only if mode = 'fixed'
%
%  Typical curvature range: 0.01–0.5 m/km². A deadband of 0.02 is
%  ~5–10% of signal.
%
%  COMMON RANGE: 0.005–0.1 m/km²
%  DEFAULT: 0.02
% -------------------------------------------------------------------------
cfg.segmentation.deadbandAmplitude_fixed = 0.02;

% -------------------------------------------------------------------------
% ★  Percentile deadband factor — used only if mode = 'percentile'
%
%  Deadband = factor × P80(|curvature|).
%  Smaller (0.05): tighter → more transitions.
%  Larger (0.2): wider → fewer transitions.
%
%  COMMON RANGE: 0.03–0.30
%  DEFAULT: 0.10
% -------------------------------------------------------------------------
cfg.segmentation.deadbandFactor_percentile = 0.10;

% -------------------------------------------------------------------------
% ★★  MAD deadband multiplier (σ units) — used only if mode = 'mad'
%
%  Deadband = deadbandSigmas × robust_σ,
%    where robust_σ = 1.4826 × MAD(curvature).
%
%  STATISTICAL CONTEXT: 1.5σ ≈ 87% confidence interval.
%
%  IF TOO SMALL (< 0.5): Almost all curvature is "significant" → dozens
%    of small, unstable domains.
%  IF TOO LARGE (> 2.5): Only extreme curvature survives → misses subtle
%    but real transitions.
%
%  ADJUST FOR:
%    Flat, low-contrast trenches (Cascadia): reduce to 1.0
%    Rough, high-contrast trenches (Izu-Bonin): increase to 2.0
%
%  COMMON RANGE: 0.5–3.0
%  DEFAULT: 1.5
% -------------------------------------------------------------------------
cfg.segmentation.deadbandSigmas_mad = 1.5;

% -------------------------------------------------------------------------
% ★★  Persistence length (km)
%
%  WHAT IT DOES:
%    A curvature state (positive or negative, outside deadband) must
%    persist for at least this distance to be confirmed as a "run".
%    Shorter excursions are treated as transient noise.
%
%  PHYSICAL MEANING: Real morphologic domains extend 10–50+ km along
%    trench axes (Clift & Vannucchi, 2004). Individual basins are
%    typically 10–20 km wide.
%
%  IF TOO SMALL (< 5 km): Noise triggers spurious state changes.
%  IF TOO LARGE (> 30 km): Individual basins become invisible.
%
%  COMMON RANGE: 3–50 km
%  DEFAULT: 10 km
% -------------------------------------------------------------------------
cfg.segmentation.persistenceLength_km = 10;

% -------------------------------------------------------------------------
% ★★  Minimum domain length after merging (km)
%
%  Domains shorter than this are merged with neighbours.
%
%  RULE: Should be ≥ 2 × persistenceLength_km.
%
%  COMMON RANGE: 5–100 km
%  DEFAULT: 20 km
% -------------------------------------------------------------------------
cfg.segmentation.minDomainLength_km = 40;

% -------------------------------------------------------------------------
% ★  Transition confidence scoring weights [sharpness, dwell, amplitude]
%
%  COMPONENTS:
%    Sharpness — How abruptly curvature changes at the boundary.
%    Dwell     — How long the state persists on each side (geometric mean).
%    Amplitude — Mean curvature contrast across the boundary.
%
%  Weights are normalised internally, so [0.4, 0.3, 0.3] and [4, 3, 3]
%  produce identical results.
%
%  DEFAULT: [0.4, 0.3, 0.3]  (slightly emphasises sharpness)
% -------------------------------------------------------------------------
cfg.segmentation.confidenceWeights = [0.5, 0.2, 0.3];

% -------------------------------------------------------------------------
% ★  Maximum number of domain boundaries
%
%  Hard upper limit to prevent over-segmentation.
%
%  COMMON RANGE: 5–100
%  DEFAULT: 30
% -------------------------------------------------------------------------
cfg.segmentation.maxBoundaries = 25;

% -------------------------------------------------------------------------
% ★  Pre-smoothing window for curvature before segmentation (km)
%
%  Gaussian smoothing of curvature before the deadband/persistence
%  analysis. Reduces pixel-scale jitter. Set to 0 to disable.
%
%  NOTE: This is separate from derivWindow_km (applied during
%  differentiation). The two smoothings are sequential.
%
%  COMMON RANGE: 0–20 km
%  DEFAULT: 5 km
% -------------------------------------------------------------------------
cfg.segmentation.smoothWindow_km = 15;

% -------------------------------------------------------------------------
% ★★  Rolling-window MAD deadband width (km)
%
%  WHAT IT DOES:
%    In deadbandMode = 'mad', the noise floor is estimated locally using a
%    rolling window of this width rather than a single global MAD.  A wider
%    window produces a smoother, less data-adaptive deadband; a narrower
%    window allows the deadband to follow local noise variations.
%
%  PHYSICAL MEANING: Roughly matches the spatial scale over which you
%    expect the curvature noise floor to be stationary (one tectonic
%    segment ≈ 50–200 km).
%
%  IF TOO SMALL (< 50 km): Deadband fluctuates rapidly → unstable domain
%    boundaries in low-curvature sections.
%  IF TOO LARGE (> 300 km): Approaches a global estimate, losing local
%    adaptivity.
%
%  COMMON RANGE: 50–200 km
%  DEFAULT: 100 km
% -------------------------------------------------------------------------
cfg.segmentation.deadbandWindow_km = 100;

% -------------------------------------------------------------------------
% ★★  Dwell normalisation scale (km)
%
%  WHAT IT DOES:
%    Controls the spatial scale used to normalise dwell length in the
%    transition confidence score.  Dwell lengths are divided by this scale
%    before entering the [0, 1] score.  Larger values reduce the weight
%    given to long-dwelling states relative to short-dwelling ones.
%
%  PHYSICAL MEANING: Should be set to the expected length of a "typical"
%    morphologic domain (one forearc basin or inter-canyon flat).
%
%  IF TOO SMALL (< 30 km): Long-dwelling states receive disproportionately
%    high confidence scores.
%  IF TOO LARGE (> 300 km): All dwell lengths appear short → confidence
%    scores cluster near zero.
%
%  COMMON RANGE: 50–200 km
%  DEFAULT: 100 km
% -------------------------------------------------------------------------
cfg.segmentation.dwellScale_km = 100;

% -------------------------------------------------------------------------
% ★★  Local amplitude normalisation radius (km)
%
%  WHAT IT DOES:
%    For the amplitude component of the confidence score, curvature
%    contrast is normalised by the local curvature amplitude within a
%    window of this half-radius around each boundary.  This makes the
%    score independent of the global curvature range.
%
%  PHYSICAL MEANING: Should be comparable to the half-width of a typical
%    domain — ~25–75 km for most trench systems.
%
%  IF TOO SMALL (< 20 km): Amplitude estimate is noisy (few points).
%  IF TOO LARGE (> 150 km): Normalisation becomes quasi-global, losing
%    local sensitivity.
%
%  COMMON RANGE: 25–100 km
%  DEFAULT: 50 km
% -------------------------------------------------------------------------
cfg.segmentation.localAmplitudeRadius_km = 50;

% -------------------------------------------------------------------------
% ★  Stability audit — parameter-sweep robustness test
%
%  After segmentation, the algorithm re-runs with ±30% perturbation of
%  {deadband, persistence, minDomainLength} over a grid (5³ = 125 combos)
%  and reports which boundaries survive across most perturbations.
%
%  Boundaries surviving < 50% of combinations are fragile and drawn as
%  dotted lines in Panel C.  If minStability > 0, fragile boundaries are
%  removed entirely.
%
%  COST: ~125 × O(n) calls, < 1 s for a 3000 km profile.
%  DEFAULT: enabled, no filtering (minStability = 0)
% -------------------------------------------------------------------------
cfg.segmentation.stabilityTest = true;          % Enable audit
cfg.segmentation.stabilityPerturbation = 0.30;  % ±30% range
cfg.segmentation.stabilityGrid = 5;             % 5 levels per parameter → 125 combos
cfg.segmentation.stabilitySnapRadius_km = 15;   % Match radius
cfg.segmentation.minStability = 0.5;            % 0 = keep all; 0.5 = keep robust only


%% ========================================================================
%  SECTION 6: SAMPLE HANDLING
%  ========================================================================
%  Column names and projection settings for the sample file.

% -------------------------------------------------------------------------
% ★  Sample file column names
%
%  Specify the exact column names for coordinates and identifiers.
%  The code also searches for common variants (case-insensitive) as
%  fallback: 'x'/'X'/'Easting', 'y'/'Y'/'Northing', etc.
%
%  COORDINATE SYSTEM: Must match the DEM. If the DEM uses projected
%  coordinates (metres), the sample file must also use metres. If
%  geographic (degrees), samples must be in degrees.
% -------------------------------------------------------------------------
cfg.samples.colX    = 'Longitude';      % X-coordinate column name
cfg.samples.colY    = 'Latitude';      % Y-coordinate column name
cfg.samples.colName = 'ID';   % Sample identifier column name
%  Additional name fallbacks searched automatically:
%    'OriginalID', 'ID', 'SampleID', 'Station'

% -------------------------------------------------------------------------
% ★  Show sample name labels on plots (logical)
%  DEFAULT: true
% -------------------------------------------------------------------------
cfg.samples.showNames = true;


%% ========================================================================
%  SECTION 7: FIGURE FORMATTING
%  ========================================================================
%  Visual parameters for the 5-panel master figure and QC diagnostic.

% -------------------------------------------------------------------------
% Figure dimensions [width, height] in pixels
%  Default aspect ratio (1.45:1) optimised for 5-row vertical layout.
%  For journal submission at 300 DPI: 1600 px = 5.3 inches.
% -------------------------------------------------------------------------
cfg.plot.figureSize = [1600, 1100];

% -------------------------------------------------------------------------
% Font size hierarchy (points)
%  title > label > tick > annotation
% -------------------------------------------------------------------------
cfg.plot.titleFontSize = 11;
cfg.plot.labelFontSize = 10;
cfg.plot.tickFontSize  = 9;
cfg.plot.annotFontSize = 8;

% -------------------------------------------------------------------------
% Sample marker size (points²)
%  50 for main figure, 40 for diagnostic (smaller panels).
% -------------------------------------------------------------------------
cfg.plot.sampleMarkerSize            = 50;
cfg.plot.sampleMarkerSize_diagnostic = 40;

% -------------------------------------------------------------------------
% Plot downsampling factor (integer ≥ 1)
%  Plot every Nth point for display performance.
%  Does NOT affect computation — only visualisation.
%  1 = plot every point, 5 = plot 20% of points.
% -------------------------------------------------------------------------
cfg.plot.downsampleFactor = 5;

% -------------------------------------------------------------------------
% Figure title prefix
%  Prepended to auto-generated figure titles.
% -------------------------------------------------------------------------
cfg.plot.titlePrefix = '';  % Empty for publication (no figure suptitle)


%% ========================================================================
%  SECTION 8: OUTPUT
%  ========================================================================

% -------------------------------------------------------------------------
% Figure export formats (cell array)
%  Supported: 'png', 'svg', 'pdf', 'eps', 'fig', 'tif'
% -------------------------------------------------------------------------
cfg.output.formats = {'png', 'svg', 'pdf'};

% -------------------------------------------------------------------------
% PNG resolution (DPI)
%  300 = publication quality, 150 = web, 600 = high-quality print.
% -------------------------------------------------------------------------
cfg.output.resolution_dpi = 600;  % Publication quality

% -------------------------------------------------------------------------
% Transparent background for PNG (logical)
%  Useful for compositing in presentations or posters.
% -------------------------------------------------------------------------
cfg.output.transparent = true;

% -------------------------------------------------------------------------
% Verbose console output (logical)
%  Prints progress messages, parameter summaries, and timings.
% -------------------------------------------------------------------------
cfg.output.verbose = true;


%% ========================================================================
%  APPENDIX: INTERNAL CONSTANTS — Do Not Change
%  ========================================================================
%  These are documented here for completeness. They are mathematical
%  constants, safety guards, or empirically validated algorithmic
%  invariants embedded in the source functions. Changing them requires
%  deep understanding of the code internals.
%
%  CONSTANT                          VALUE     FILE : LINE      JUSTIFICATION
%  ─────────────────────────────────────────────────────────────────────────
%  MAD-to-σ conversion factor        1.4826    segmentDomains     Exact: 1/Φ⁻¹(3/4) for Gaussian.
%  Deadband floor                    1e-6      segmentDomains     Machine-precision guard.
%  Gaussian kernel truncation        4σ        computeTrenchBands Standard DSP; captures 99.994%.
%  Segmentation smooth multiplier    3×        computeTrenchBands Segmentation κ uses 3× wider window.
%  Tectonic derivative window        3 samples computeTrenchBands Tectonic bands already smooth.
%  Min points for analysis           100       computeTrenchBands Below 100, decomposition meaningless.
%  Min points in metric window       10        computeSinkMetrics Statistical minimum for RMS/MAD.
%  Min points for Welch PSD          20        computeSinkMetrics Welch's method minimum.
%  Min frequency bins for β fit      4         computeSinkMetrics Power-law fit minimum.
%  Min drainage basin area           10 cells  extractMainDivide  Below 10 cells = grid artefact.
%  Duplicate point removal           0.1×cell  orderPoints        10% of cell = same physical point.
%  Percentile-mode reference         P80       segmentDomains     Scales to 80th %ile of |curvature|.
%  Max dwell fraction                0.5       segmentDomains     Max dwell = half profile length.
%  Sharpness divisor                 10        segmentDomains     Heuristic normalisation.
%  TI reference lines                0.5, 0.7  plotTrenchBands    Empirical "moderate/high" thresholds.
%  Sample label halo                 12/11pt   plotTrenchBands    White background + black foreground.
%  Y-limit percentile × margin       various   plotTrenchBands    Adaptive axis limits.
%  ─────────────────────────────────────────────────────────────────────────

end
