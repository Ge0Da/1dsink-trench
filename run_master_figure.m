%% run_master_figure.m
% =========================================================================
%  1D Along-Strike Sink Master Figure
% =========================================================================
%
%  This script produces the 1D Along-Strike Sink Master Figure (5 panels) 
%  and the QC Diagnostic figure (6 panels) from a bathymetric DEM and sample
%  information.
%
%  PRODUCTION CHAIN:
%    Step 1:  Load configuration (config_master_figure.m)
%    Step 2:  Validate inputs (DEM + samples are MANDATORY)
%    Step 3:  Load DEM
%    Step 4:  Extract trench axis from DEM
%    Step 5:  Load sample table
%    Step 6:  Project samples onto trench axis
%    Step 7:  Multi-band spectral decomposition
%    Step 8:  Produce 1DSink Master Figure (5-panel)
%    Step 9:  Produce QC Diagnostic Figure (6-panel)
%    Step 10: Save figures
%    Step 11: Save data products (MAT, CSV)
%
%  USAGE:
%    1. Edit config_master_figure.m (file paths, parameters)
%    2. Run this script: >> run_master_figure
%
%  PREREQUISITES:
%    - TopoToolbox on MATLAB path
%    - All utility functions in the utils/ subdirectory (see README.md)
%    - DEM file readable by GRIDobj (GeoTIFF, NetCDF, or MAT)
%    - Sample file readable by readtable (CSV, XLSX)
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%  Author: David Fernández-Blanco
%  Date:   2026
%  Version: 1.1
% =========================================================================


%% ========================================================================
%  STEP 0: ENVIRONMENT — add source directories to MATLAB path
%  ========================================================================
%  Adds this directory and all subdirectories (including utils/) to the
%  MATLAB path. This makes the package fully self-contained — no external
%  project dependencies.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

addpath(genpath(scriptDir));

% --- Verify critical functions are reachable ---
criticalFunctions = {'resampleToCoarse', 'computeFlowRouting', ...
    'extractMainDivide', 'computeTrenchBands', 'segmentDomains', ...
    'extractTrenchAxis', 'projectSamplesOnAxis', 'trimBands'};
missing = criticalFunctions(cellfun(@(f) exist(f,'file')==0, criticalFunctions));
if ~isempty(missing)
    error('run_master_figure:missingFunctions', ...
        ['Cannot find required functions: %s\n\n' ...
         'All utility functions must be in the utils/ subdirectory.\n' ...
         'See README.md for the complete file list.'], ...
        strjoin(missing, ', '));
end


%% ========================================================================
%  STEP 1: LOAD CONFIGURATION
%  ========================================================================
%  All parameters come from config_master_figure.m.
%  Edit that file — NOT this driver — to change behaviour.

fprintf('\n');
fprintf('================================================================\n');
fprintf('  1D SINK MASTER FIGURE — Production Driver\n');
fprintf('================================================================\n\n');

cfg = config_master_figure();

verbose = cfg.output.verbose;
if verbose, fprintf('Step 1: Configuration loaded.\n'); end


%% ========================================================================
%  STEP 2: VALIDATE INPUTS — files must exist
%  ========================================================================

% Resolve paths relative to the driver location
demPath    = resolveRelativePath(cfg.demFile, scriptDir);
samplePath = resolveRelativePath(cfg.sampleFile, scriptDir);
outputDir  = resolveRelativePath(cfg.outputDir, scriptDir);

% DEM is mandatory
if ~isfile(demPath)
    error('run_master_figure:demNotFound', ...
        'DEM file not found:\n  %s\nEdit cfg.demFile in config_master_figure.m', demPath);
end

% Samples are MANDATORY (RE-1: not optional)
if ~isfile(samplePath)
    error('run_master_figure:samplesNotFound', ...
        ['Sample file not found:\n  %s\n' ...
         'Samples are REQUIRED for this production figure.\n' ...
         'Edit cfg.sampleFile in config_master_figure.m'], samplePath);
end

% Create output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    if verbose, fprintf('  Created output directory: %s\n', outputDir); end
end

if verbose
    fprintf('Step 2: Input validation passed.\n');
    fprintf('        DEM:     %s\n', demPath);
    fprintf('        Samples: %s\n', samplePath);
    fprintf('        Output:  %s\n', outputDir);
end


%% ========================================================================
%  STEP 3: LOAD DEM
%  ========================================================================
if verbose, fprintf('\nStep 3: Loading DEM...\n'); end
tStart = tic;

DEM = GRIDobj(demPath);

if verbose
    fprintf('  Grid size:     %d × %d\n', DEM.size(1), DEM.size(2));
    fprintf('  Cell size:     %.1f m\n', DEM.cellsize);
    fprintf('  Depth range:   %.0f to %.0f m\n', ...
        min(DEM.Z(:), [], 'omitnan'), max(DEM.Z(:), [], 'omitnan'));
    fprintf('  Loaded in %.1f s\n', toc(tStart));
end


%% ========================================================================
%  STEP 4: EXTRACT TRENCH AXIS
%  ========================================================================
%  TRIM-INVARIANCE: Always extract the FULL axis first.  Segmentation and
%  spectral decomposition run on the full profile so that results are
%  independent of the display window.  Trimming happens AFTER analysis,
%  before plotting.
if verbose, fprintf('\nStep 4: Extracting trench axis...\n'); end
tStart = tic;

trenchAxis = extractTrenchAxis(DEM, ...
    'coarseRes_m',              cfg.axis.coarseRes_m, ...
    'minAreaKm2',               cfg.axis.minAreaKm2, ...
    'localRefineRadius_cells',  cfg.axis.localRefineRadius_cells, ...
    'smoothWindow_samples',     cfg.axis.smoothWindow_samples, ...
    'axisTrim_km',              [], ...      % <-- NO TRIM (full profile)
    'refineToFullRes',          cfg.axis.refineToFullRes, ...
    'verbose',                  verbose);

if verbose, fprintf('  Axis extraction: %.1f s\n', toc(tStart)); end


%% ========================================================================
%  STEP 5: LOAD SAMPLES
%  ========================================================================
if verbose, fprintf('\nStep 5: Loading sample table...\n'); end

sampleTable = readtable(samplePath);
if verbose
    fprintf('  %d samples, %d columns\n', height(sampleTable), width(sampleTable));
    fprintf('  Columns: %s\n', strjoin(sampleTable.Properties.VariableNames, ', '));
end

% --- Reproject sample coordinates if needed ---
% Sample files store Longitude/Latitude (degrees); the DEM and trench axis
% are in a projected CRS (metres).  Detect the mismatch and convert.
% This mirrors what ModuleA_Main_Driver does at line 494-501 after A01.
sp_x = sampleTable.(cfg.samples.colX);
sp_y = sampleTable.(cfg.samples.colY);

if all(abs(sp_x) <= 360) && all(abs(sp_y) <= 90) ...
        && max(abs(trenchAxis.x(:))) > 10000
    if verbose
        fprintf('  Samples in degrees, axis in metres — reprojecting...\n');
    end
    [sp_x, sp_y, crsMethod] = reprojectSamplesToGrid(sp_x, sp_y, demPath);
    sampleTable.(cfg.samples.colX) = sp_x;
    sampleTable.(cfg.samples.colY) = sp_y;
    if verbose
        fprintf('  Done (%s). Range: x %.0f–%.0f, y %.0f–%.0f\n', ...
            crsMethod, min(sp_x), max(sp_x), min(sp_y), max(sp_y));
    end
end


%% ========================================================================
%  STEP 6: PROJECT SAMPLES ONTO TRENCH AXIS
%  ========================================================================
if verbose, fprintf('\nStep 6: Projecting samples onto axis...\n'); end

projectedSamples = projectSamplesOnAxis(sampleTable, trenchAxis, ...
    'colX',      cfg.samples.colX, ...
    'colY',      cfg.samples.colY, ...
    'colName',   cfg.samples.colName, ...
    'showNames', cfg.samples.showNames, ...
    'verbose',   verbose);


%% ========================================================================
%  STEP 7: MULTI-BAND SPECTRAL DECOMPOSITION (full profile)
%  ========================================================================
if verbose, fprintf('\nStep 7: Computing multi-band decomposition...\n'); end
tStart = tic;

bands = computeTrenchBands(trenchAxis.s_m, trenchAxis.z, ...
    'resampleSpacing_m',  cfg.signal.resampleSpacing_m, ...
    'interpMethod',       cfg.signal.interpMethod, ...
    'edgeTrim_km',        cfg.signal.edgeTrim_km, ...
    'derivWindow_km',     cfg.signal.derivWindow_km, ...
    'lambda_km',          cfg.signal.lambda_km, ...
    'useDepth',           true);

if verbose
    fprintf('  Uniform grid: %d points at %.0f m spacing\n', ...
        bands.nPoints, cfg.signal.resampleSpacing_m);
    fprintf('  Profile length: %.1f km\n', bands.length_km);
    fprintf('  Decomposition: %.1f s\n', toc(tStart));
end


%% ========================================================================
%  STEP 7b: SEGMENTATION ON FULL PROFILE
%  ========================================================================
%  segmentDomains v2.0 is spatially stationary: results are independent
%  of the profile extent.  We still run on the full profile so that
%  edge-adjacent domains get the most context available.
if verbose, fprintf('\nStep 7b: Running segmentation on full profile...\n'); end
tStart = tic;

[segDomains, segBounds, segDiag] = segmentDomains( ...
    bands.s_km, bands.d2D_2_20_mpkm2, bands.edgeMask, ...
    'minDomainLength_km',   cfg.segmentation.minDomainLength_km, ...
    'persistenceLength_km', cfg.segmentation.persistenceLength_km, ...
    'deadbandMode',         cfg.segmentation.deadbandMode, ...
    'deadbandAmplitude',    cfg.segmentation.deadbandAmplitude_fixed, ...
    'deadbandFactor',       cfg.segmentation.deadbandFactor_percentile, ...
    'deadbandSigmas',       cfg.segmentation.deadbandSigmas_mad, ...
    'deadbandWindow_km',    cfg.segmentation.deadbandWindow_km, ...
    'confidenceWeights',    cfg.segmentation.confidenceWeights, ...
    'dwellScale_km',        cfg.segmentation.dwellScale_km, ...
    'localAmplitudeRadius_km', cfg.segmentation.localAmplitudeRadius_km, ...
    'maxBoundaries',        cfg.segmentation.maxBoundaries, ...
    'smoothWindow_km',      cfg.segmentation.smoothWindow_km, ...
    'stabilityTest',        cfg.segmentation.stabilityTest, ...
    'stabilityPerturbation',cfg.segmentation.stabilityPerturbation, ...
    'stabilityGrid',        cfg.segmentation.stabilityGrid, ...
    'stabilitySnapRadius_km',cfg.segmentation.stabilitySnapRadius_km, ...
    'minStability',         cfg.segmentation.minStability, ...
    'verbose',              verbose);

precomputedSeg = struct('domains', segDomains, ...
                        'boundaries', segBounds, ...
                        'diagnostics', segDiag);

if verbose
    fprintf('  %d domains, %d boundaries (%.1f s)\n', ...
        numel(segDomains), numel(segBounds), toc(tStart));
end


%% ========================================================================
%  STEP 7b-conv: CONVERGENCE SCORING ON FULL PROFILE
%  ========================================================================
%  Compute continuous multi-proxy convergence scores for each boundary.
%  Uses the geometric mean of contrast ratios from curvature, slope, and
%  roughness — no thresholding, no discrete scores.
if verbose, fprintf('\nStep 7b-conv: Computing convergence scores...\n'); end

convWindowHalfWidth_km = 20;   % half-width for before/after comparison
[convComposites, convDetails] = convergenceScore( ...
    segBounds, bands.s_km, ...
    bands.d2D_2_20_mpkm2, ...   % curvature signal
    bands.dD_2_10_mpkm, ...     % slope signal
    bands.D_1_10_m, ...         % roughness signal
    convWindowHalfWidth_km);

precomputedSeg.convergence = convDetails;

if verbose
    fprintf('  %d boundaries scored, composite CR range [%.3f, %.3f]\n', ...
        numel(convComposites), min(convComposites), max(convComposites));
end


%% ========================================================================
%  STEP 7c: TRIM BANDS FOR DISPLAY (if requested)
%  ========================================================================
%  The axisTrim_km parameter now controls the DISPLAY window only.
%  Segmentation results from Step 7b are passed through unchanged and
%  clipped to the display range inside plotTrenchBands_1DSink.
displayTrim = cfg.axis.axisTrim_km;
if ~isempty(displayTrim)
    if verbose
        fprintf('\nStep 7c: Trimming bands for display [%.0f, %.0f] km...\n', ...
            displayTrim(1), displayTrim(2));
    end
    bands = trimBands(bands, displayTrim);
    if verbose
        fprintf('  Display range: %.1f km (%d points)\n', ...
            bands.length_km, bands.nPoints);
    end
end


%% ========================================================================
%  STEP 8: PRODUCE 1D SINK MASTER FIGURE (5-panel)
%  ========================================================================
if verbose, fprintf('\nStep 8: Generating 1D Sink Master Figure...\n'); end
tStart = tic;

[figMaster, masterResults] = plotTrenchBands_1DSink(bands, ...
    'title',                cfg.plot.titlePrefix, ...
    'samples',              projectedSamples, ...
    'showSampleNames',      cfg.samples.showNames, ...
    'sampleMarkerSize',     cfg.plot.sampleMarkerSize, ...
    'downsamplePlot',       cfg.plot.downsampleFactor, ...
    ... % --- Sink metrics ---
    'windowSize_km',        cfg.metrics.windowSize_km, ...
    'windowOverlap',        cfg.metrics.windowOverlap, ...
    'w_slope',              cfg.metrics.w_slope, ...
    'w_rough',              cfg.metrics.w_rough, ...
    'w_curv',               cfg.metrics.w_curv, ...
    'normMode',             cfg.metrics.normMode, ...
    'slopeLo_mpkm',         cfg.metrics.slopeLo_mpkm, ...
    'slopeHi_mpkm',         cfg.metrics.slopeHi_mpkm, ...
    'roughLo_m',            cfg.metrics.roughLo_m, ...
    'roughHi_m',            cfg.metrics.roughHi_m, ...
    'curvLo_mpkm2',         cfg.metrics.curvLo_mpkm2, ...
    'curvHi_mpkm2',         cfg.metrics.curvHi_mpkm2, ...
    'pctileLow',            cfg.metrics.pctileLow, ...
    'pctileHigh',           cfg.metrics.pctileHigh, ...
    'betaDetrend',          cfg.metrics.betaDetrend, ...
    'betaFitBand',          cfg.metrics.betaFitBand_cyc_km, ...
    ... % --- Domain segmentation ---
    'minDomainLength_km',   cfg.segmentation.minDomainLength_km, ...
    'persistenceLength_km', cfg.segmentation.persistenceLength_km, ...
    'deadbandMode',         cfg.segmentation.deadbandMode, ...
    'deadbandAmplitude',    cfg.segmentation.deadbandAmplitude_fixed, ...
    'deadbandFactor',       cfg.segmentation.deadbandFactor_percentile, ...
    'deadbandSigmas',       cfg.segmentation.deadbandSigmas_mad, ...
    'deadbandWindow_km',    cfg.segmentation.deadbandWindow_km, ...
    'confidenceWeights',    cfg.segmentation.confidenceWeights, ...
    'dwellScale_km',        cfg.segmentation.dwellScale_km, ...
    'localAmplitudeRadius_km', cfg.segmentation.localAmplitudeRadius_km, ...
    'maxBoundaries',        cfg.segmentation.maxBoundaries, ...
    'segmentSmoothWindow_km', cfg.segmentation.smoothWindow_km, ...
    ... % --- Stability audit ---
    'stabilityTest',          cfg.segmentation.stabilityTest, ...
    'stabilityPerturbation',  cfg.segmentation.stabilityPerturbation, ...
    'stabilityGrid',          cfg.segmentation.stabilityGrid, ...
    'stabilitySnapRadius_km', cfg.segmentation.stabilitySnapRadius_km, ...
    'minStability',           cfg.segmentation.minStability, ...
    ... % --- Pre-computed segmentation (trim-invariant) ---
    'precomputedSegmentation', precomputedSeg);

% Apply figure formatting from config
set(figMaster, 'Position', [30 30 cfg.plot.figureSize(1) cfg.plot.figureSize(2)]);

if verbose, fprintf('  Master figure: %.1f s\n', toc(tStart)); end


%% ========================================================================
%  STEP 9: PRODUCE QC DIAGNOSTIC FIGURE (6-panel)
%  ========================================================================
if verbose, fprintf('\nStep 9: Generating QC Diagnostic Figure...\n'); end
tStart = tic;

[figDiag, diagMetrics] = plotSinkMetricsDiagnostic(bands, ...
    'title',            'Sink Metrics Diagnostic — RMS, MAD, β Computation', ...
    'windowSize_km',    cfg.metrics.windowSize_km, ...
    'windowOverlap',    cfg.metrics.windowOverlap, ...
    'betaDetrend',      cfg.metrics.betaDetrend, ...
    'betaFitBand',      cfg.metrics.betaFitBand_cyc_km, ...
    'w_slope',          cfg.metrics.w_slope, ...
    'w_rough',          cfg.metrics.w_rough, ...
    'w_curv',           cfg.metrics.w_curv, ...
    'normMode',         cfg.metrics.normMode, ...
    'slopeLo_mpkm',     cfg.metrics.slopeLo_mpkm, ...
    'slopeHi_mpkm',     cfg.metrics.slopeHi_mpkm, ...
    'roughLo_m',        cfg.metrics.roughLo_m, ...
    'roughHi_m',        cfg.metrics.roughHi_m, ...
    'curvLo_mpkm2',     cfg.metrics.curvLo_mpkm2, ...
    'curvHi_mpkm2',     cfg.metrics.curvHi_mpkm2, ...
    'pctileLow',        cfg.metrics.pctileLow, ...
    'pctileHigh',       cfg.metrics.pctileHigh, ...
    'samples',          projectedSamples, ...
    'showSampleNames',  cfg.samples.showNames, ...
    'sampleMarkerSize', cfg.plot.sampleMarkerSize_diagnostic);

if verbose, fprintf('  Diagnostic figure: %.1f s\n', toc(tStart)); end


%% ========================================================================
%  STEP 10: SAVE FIGURES
%  ========================================================================
if verbose, fprintf('\nStep 10: Saving figures...\n'); end

saveFigure(figMaster, outputDir, '1DSink_MasterFigure', cfg.output.formats, ...
    'resolution', cfg.output.resolution_dpi, ...
    'transparent', cfg.output.transparent, ...
    'verbose', verbose);

saveFigure(figDiag, outputDir, '1DSink_Diagnostic', cfg.output.formats, ...
    'resolution', cfg.output.resolution_dpi, ...
    'transparent', cfg.output.transparent, ...
    'verbose', verbose);


%% ========================================================================
%  STEP 11: SAVE DATA PRODUCTS
%  ========================================================================
if verbose, fprintf('\nStep 11: Saving data products...\n'); end

% ---- MAT file with everything ----
matFile = fullfile(outputDir, '1DSink_MasterFigure_results.mat');
save(matFile, 'trenchAxis', 'projectedSamples', 'bands', ...
    'masterResults', 'diagMetrics', 'precomputedSeg', ...
    'convComposites', 'convDetails', 'cfg', '-v7.3');
if verbose, fprintf('  Saved: %s\n', matFile); end

% ---- CSV: axis coordinates + band depths ----
nPts = bands.nPoints;
bandCSV = table( ...
    bands.s_km(:), bands.depth_raw_m(:), ...
    bands.D_lt10_m(:), bands.D_2_10_m(:), bands.D_5_10_m(:), ...
    bands.D_10_50_m(:), bands.D_50_200_m(:), bands.D_200_1000_m(:), bands.D_gt1000_m(:), ...
    'VariableNames', { ...
        'Distance_km', 'Depth_raw_m', ...
        'D_lt10km_m', 'D_2_10km_m', 'D_5_10km_m', ...
        'D_10_50km_m', 'D_50_200km_m', 'D_200_1000km_m', 'D_gt1000km_m'});
csvFile = fullfile(outputDir, '1DSink_band_decomposition.csv');
writetable(bandCSV, csvFile);
if verbose, fprintf('  Saved: %s\n', csvFile); end

% ---- CSV: projected sample positions ----
sampleCSV = table( ...
    projectedSamples.name, ...
    projectedSamples.x, projectedSamples.y, ...
    projectedSamples.s_km, projectedSamples.d_km, projectedSamples.dist_km, ...
    'VariableNames', {'Name', 'X', 'Y', 'AlongAxis_km', 'LateralOffset_km', 'DistToAxis_km'});
sampleCSVFile = fullfile(outputDir, '1DSink_projected_samples.csv');
writetable(sampleCSV, sampleCSVFile);
if verbose, fprintf('  Saved: %s\n', sampleCSVFile); end


%% ========================================================================
%  COMPLETE
%  ========================================================================
fprintf('\n================================================================\n');
fprintf('  1D SINK MASTER FIGURE — COMPLETE\n');
fprintf('================================================================\n');
fprintf('  Trench axis:    %.1f km (%d nodes)\n', trenchAxis.length_km, trenchAxis.nNodes);
fprintf('  Samples:        %d projected\n', projectedSamples.nSamples);
fprintf('  Bands (full):   %d points at %.0f m spacing\n', bands.nPoints, cfg.signal.resampleSpacing_m);
if ~isempty(displayTrim)
    fprintf('  Display trim:   [%.0f, %.0f] km (%.1f km shown)\n', ...
        displayTrim(1), displayTrim(2), bands.length_km);
end
fprintf('  Segmentation:   %d domains, %d boundaries\n', ...
    numel(segDomains), numel(segBounds));
if ~isempty(masterResults) && isfield(masterResults, 'domains')
    fprintf('  Domains (view): %d visible\n', numel(masterResults.domains));
end
fprintf('  Output:         %s\n', outputDir);
fprintf('================================================================\n\n');


%% ========================================================================
%  LOCAL HELPER: resolve relative paths
%  ========================================================================
function absPath = resolveRelativePath(relPath, baseDir)
%RESOLVERELATIVEPATH  Convert relative path to absolute, using baseDir.
%
%   If relPath is already absolute, return it unchanged.
%   Otherwise, prepend baseDir.
    relPath = char(relPath);
    if isempty(relPath)
        absPath = baseDir;
        return;
    end
    % Check if already absolute (starts with / on Unix, X:\ on Windows)
    if relPath(1) == '/' || relPath(1) == '~' || ...
            (numel(relPath) >= 2 && relPath(2) == ':')
        absPath = relPath;
    else
        absPath = fullfile(baseDir, relPath);
    end
end
