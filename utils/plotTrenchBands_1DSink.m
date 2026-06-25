function [fig, results] = plotTrenchBands_1DSink(bands, opts)
%plotTrenchBands_1DSink Master 1DSink Routing Analysis Figure (5 profiles)
%
%   fig = plotTrenchBands_1DSink(bands)
%   [fig, results] = plotTrenchBands_1DSink(bands, 'Name', value, ...)
%
%   Creates a comprehensive 5-profile figure for Source-to-Sink analysis:
%
%   PROFILE A - SINK GEOMETRY (Architecture):
%     Left axis: LP>10 (thick blue depth envelope) + samples
%     Right axis: D_{2-10} as filled band centered at 0
%     Purpose: Show depth envelope + residual relief
%
%   PROFILE B - ROUTING GRADIENTS:
%     Thin line: Slope of LP>10 trend (envelope gradient)
%     Thick line: Moving median(|slope|) of D_{2-10} (local variability)
%     Purpose: Short-scale routing vs long-wavelength redistribution potential
%
%   PROFILE C - SINK SEGMENTATION (Curvature Domains):
%     κ_{5-10} smoothed (domain curvature)
%     Boundaries from κ_{2-20} (vertical lines)
%     Domain fill shading by sign (concave-up vs concave-down)
%     Annotation: mean roughness per domain
%     Purpose: Stable domain interpretation for ponding vs bypass
%
%   PROFILE D - SINK STRENGTH (Roughness Metrics):
%     3-row heatmap: RMS, MAD, spectral slope β
%     RMS/MAD computed on D_{2-10} (band-passed)
%     β computed on DETRENDED raw depth (raw - LP_50) for physical validity
%     Purpose: Where pathway is "shredding-prone" at 2-10 km scales
%
%   PROFILE E - TRAPPING INDEX (Composite):
%     TI = weighted(low_slope + high_roughness + concave_curvature)
%     Range 0-1: higher = higher trapping potential
%     Purpose: Synthesize sink potential quantitatively
%
%   VERSION 4.0 - Complete 1DSink Master Figure redesign
%
%   INPUTS:
%       bands   - Output from computeTrenchBands()
%
%   NAME-VALUE OPTIONS:
%       'title'             - Figure title (default: '1DSink Master Analysis')
%       'samples'           - Sample table/struct with s_axis, ID fields
%       'showSampleNames'   - Show sample labels (default: true)
%       'sampleMarkerSize'  - Marker size (default: 50)
%       'windowSize_km'     - Moving window for metrics (default: 50)
%       'minDomainLength_km' - Min domain size (default: 20)
%       'persistenceLength_km' - Min run-length for state (default: 10)
%       'deadbandMode'      - 'mad', 'percentile', or 'fixed' (default: 'mad')
%       'downsamplePlot'    - Downsample factor (default: 5)
%       'outputPath'        - Save PNG/PDF if provided
%       'w_slope'           - TI weight for slope (default: 1.0)
%       'w_rough'           - TI weight for roughness (default: 1.2)
%       'w_curv'            - TI weight for curvature (default: 0.8)
%
%   OUTPUTS:
%       fig     - Figure handle
%       results - Struct with computed metrics, domains, boundaries
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 4.0

arguments
    bands struct
    opts.title char = ''
    opts.samples = []
    opts.showSampleNames (1,1) logical = true
    opts.sampleMarkerSize (1,1) double = 50
    opts.panels (1,:) cell = {'a','b','c','d','e'}   % subset of profiles to render
    opts.windowSize_km (1,1) double = 50
    opts.minDomainLength_km (1,1) double = 20
    opts.persistenceLength_km (1,1) double = 10
    opts.deadbandMode (1,:) char {mustBeMember(opts.deadbandMode, {'mad', 'percentile', 'fixed'})} = 'mad'
    opts.downsamplePlot (1,1) double = 5
    opts.outputPath (1,:) char = ''
    opts.w_slope (1,1) double = 1.0
    opts.w_rough (1,1) double = 1.2
    opts.w_curv (1,1) double = 0.8
    % --- Tier B passthrough: computeSinkMetrics ---
    opts.windowOverlap (1,1) double = 0.8
    opts.pctileLow (1,1) double = 5
    opts.pctileHigh (1,1) double = 95
    opts.betaDetrend (1,:) char {mustBeMember(opts.betaDetrend, {'lp50', 'lp1000'})} = 'lp50'
    opts.betaFitBand (1,2) double = [0.02, 0.5]
    % --- Tier B passthrough: segmentDomains ---
    opts.deadbandAmplitude (1,1) double = 0.02
    opts.deadbandFactor (1,1) double = 0.10
    opts.deadbandSigmas (1,1) double = 1.5
    opts.confidenceWeights (1,3) double = [0.4, 0.3, 0.3]
    opts.maxBoundaries (1,1) double = 30
    opts.segmentSmoothWindow_km (1,1) double = 5
    % --- Tier B passthrough: stability audit ---
    opts.stabilityTest (1,1) logical = true
    opts.stabilityPerturbation (1,1) double = 0.30
    opts.stabilityGrid (1,1) double = 5
    opts.stabilitySnapRadius_km (1,1) double = 15
    opts.minStability (1,1) double = 0
    % --- Pre-computed segmentation (for trim-invariance) ---
    %   When provided, plotTrenchBands_1DSink skips segmentDomains and uses
    %   these results instead (clipped to the current display range).
    %   This ensures that domain boundaries are identical regardless of
    %   whether the display range is the full profile or a trimmed subset.
    %   Struct with fields: .domains, .boundaries, .diagnostics
    opts.precomputedSegmentation struct = struct()
    % --- v1.0 parameters (passed through from run_master_figure) ---
    opts.normMode (1,:) char = 'absolute'
    opts.slopeLo_mpkm (1,1) double = 0
    opts.slopeHi_mpkm (1,1) double = 50
    opts.roughLo_m (1,1) double = 0
    opts.roughHi_m (1,1) double = 30
    opts.curvLo_mpkm2 (1,1) double = 0
    opts.curvHi_mpkm2 (1,1) double = 0.5
    opts.deadbandWindow_km (1,1) double = 100
    opts.dwellScale_km (1,1) double = 100
    opts.localAmplitudeRadius_km (1,1) double = 50
end

%% Validate bands structure
hasV3 = isfield(bands, 'D_2_10_m') && isfield(bands, 'd2D_5_10_mpkm2') && isfield(bands, 'd2D_2_20_mpkm2');
if ~hasV3
    error('plotTrenchBands_1DSink:missingFields', ...
        'bands structure missing v3.0 fields. Run with updated computeTrenchBands.m');
end

%% Extract core data
s_km = bands.s_km;
ds_km = bands.ds_km;
xLim = [s_km(1), s_km(end)];        % Use actual range (not assumed 0-based)
edgeMask = bands.edgeMask;
ds = opts.downsamplePlot;

%% Compute sink metrics
metrics = computeSinkMetrics(bands, ...
    'windowSize_km', opts.windowSize_km, ...
    'windowOverlap', opts.windowOverlap, ...
    'w_slope', opts.w_slope, ...
    'w_rough', opts.w_rough, ...
    'w_curv', opts.w_curv, ...
    'pctileLow', opts.pctileLow, ...
    'pctileHigh', opts.pctileHigh, ...
    'betaDetrend', opts.betaDetrend, ...
    'betaFitBand', opts.betaFitBand);

%% Auto-compute valid display range from the data
% Valid = edgeMask (filter signals) AND ~isnan(TI_full) (windowed metrics).
% The union of both constraints gives the range where ALL panels have
% meaningful data — no hardcoded km values, no manual tuning needed.
validDisplay = edgeMask(:) & ~isnan(metrics.TI_full(:));
validIdx_disp = find(validDisplay);
if ~isempty(validIdx_disp)
    xLim = [s_km(validIdx_disp(1)), s_km(validIdx_disp(end))];
else
    xLim = [s_km(1), s_km(end)];   % fallback: show full range
end

%% Segment boundaries — use pre-computed results or compute fresh
%
%  TRIM-INVARIANCE PRINCIPLE:
%  Segmentation must produce the same boundaries regardless of display
%  range.  When pre-computed results from a full-profile run are available,
%  we clip them to the current display range instead of re-running the
%  algorithm on a trimmed window (which would change the deadband, confidence
%  normalization, and edge constraints).

hasPre = isfield(opts.precomputedSegmentation, 'domains') && ...
         isfield(opts.precomputedSegmentation, 'boundaries');

if hasPre
    % --- Use pre-computed segmentation (clip to display range) ---
    allDomains = opts.precomputedSegmentation.domains;
    allBounds  = opts.precomputedSegmentation.boundaries;
    segDiag    = opts.precomputedSegmentation.diagnostics;

    % Clip boundaries to display range
    segmentBounds = allBounds(allBounds >= xLim(1) & allBounds <= xLim(2));

    % Clip domains: keep any domain that overlaps the display range
    domains = allDomains([]);   % empty struct with same fields
    for d = 1:numel(allDomains)
        dom = allDomains(d);
        % Does this domain overlap the display window?
        if dom.s_end >= xLim(1) && dom.s_start <= xLim(2)
            % Clip domain extent to display window
            dom.s_start = max(dom.s_start, xLim(1));
            dom.s_end   = min(dom.s_end,   xLim(2));
            dom.length_km = dom.s_end - dom.s_start;
            domains(end+1) = dom; %#ok<AGROW>
        end
    end
    % Re-number domain IDs for display
    for d = 1:numel(domains), domains(d).id = d; end

    % Clip boundary-indexed vectors to match clipped boundaries
    keepBound = allBounds >= xLim(1) & allBounds <= xLim(2);
    if isfield(segDiag, 'boundaryStability') && ~isempty(segDiag.boundaryStability)
        segDiag.boundaryStability = segDiag.boundaryStability(keepBound);
    end

    % Clip convergence data if present
    hasConv = isfield(opts.precomputedSegmentation, 'convergence');
    if hasConv
        convData = opts.precomputedSegmentation.convergence;
        fnames = fieldnames(convData);
        for fi = 1:numel(fnames)
            fn = fnames{fi};
            val = convData.(fn);
            % Clip row-vectors that match the original boundary count
            if isnumeric(val) && size(val,1) == numel(allBounds)
                convData.(fn) = val(keepBound, :);
            end
        end
    else
        convData = [];
    end
else
    % --- No pre-computed results: run segmentDomains on current window ---
    %   WARNING: results will change if the display range changes.
    [domains, segmentBounds, segDiag] = segmentDomains(s_km, bands.d2D_2_20_mpkm2, edgeMask, ...
        'minDomainLength_km', opts.minDomainLength_km, ...
        'persistenceLength_km', opts.persistenceLength_km, ...
        'deadbandMode', opts.deadbandMode, ...
        'deadbandAmplitude', opts.deadbandAmplitude, ...
        'deadbandFactor', opts.deadbandFactor, ...
        'deadbandSigmas', opts.deadbandSigmas, ...
        'confidenceWeights', opts.confidenceWeights, ...
        'maxBoundaries', opts.maxBoundaries, ...
        'smoothWindow_km', opts.segmentSmoothWindow_km, ...
        'stabilityTest', opts.stabilityTest, ...
        'stabilityPerturbation', opts.stabilityPerturbation, ...
        'stabilityGrid', opts.stabilityGrid, ...
        'stabilitySnapRadius_km', opts.stabilitySnapRadius_km, ...
        'minStability', opts.minStability, ...
        'verbose', false);
    hasConv  = false;
    convData = [];
end

%% Process samples
[hasSamples, sample_s_km, sample_names, nSamp] = processSamples(opts.samples);

%% Compute additional derivatives needed
% Slope of LP>10 (envelope gradient)
LP10_slope = gradient(bands.lowpass.LP_10, ds_km);  % m/km

%% Panel selection
do_a = any(strcmp(opts.panels, 'a'));
do_b = any(strcmp(opts.panels, 'b'));
do_c = any(strcmp(opts.panels, 'c'));
do_d = any(strcmp(opts.panels, 'd'));
do_e = any(strcmp(opts.panels, 'e'));
nPanels = do_a + do_b + do_c + do_d + do_e;
if nPanels == 0
    error('plotTrenchBands_1DSink:noPanels', 'opts.panels must request at least one of a-e.');
end

%% Create figure
% Height scales with the number of requested panels so each tile keeps a
% similar aspect ratio regardless of how many profiles are drawn.
figH = max(220, round(1100 * nPanels / 5));
fig = figure('Name', '1DSink Master Analysis', ...
    'Position', [30 30 1600 figH], 'Color', 'w');

tiledlayout(nPanels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel handles (populated only for requested profiles).
ax1 = []; ax2 = []; ax3 = []; ax4 = []; ax5 = [];

%% Colors and styles
colorLP = [0.1 0.3 0.7];           % Blue for LP>10
colorResidual = [0.2 0.6 0.3];     % Green for D_2_10
colorSlopeBroad = [0.4 0.4 0.4];   % Gray for LP slope
colorSlopeLocal = [0.8 0.3 0.2];   % Red for local slope
colorCurvature = [0.5 0.2 0.6];    % Purple for curvature
colorTI = [0.9 0.5 0.1];           % Orange for TI
colorSample = [0.95 0.6 0.1];      % Orange for samples
colorBoundary = [0.7 0.2 0.2];     % Red for boundaries
colorDomainPos = [0.78 0.88 0.98];  % Light blue for bypass (concave-down)
colorDomainNeg = [0.99 0.90 0.72];  % Light amber for trap (concave-up)

% Font sizes
titleFS = 11;
labelFS = 10;
tickFS = 9;
annotFS = 8;

%% ========================================================================
%  PROFILE A: SINK GEOMETRY (Architecture)
%  ========================================================================
if do_a
ax1 = nexttile;
hold on;

% --- Step 1: RIGHT AXIS (background) — D_{2-10} residual band ---
yyaxis right;

D_residual = bands.D_2_10_m;
s_ds = s_km(1:ds:end);
D_ds = D_residual(1:ds:end);

fill([s_ds; flipud(s_ds)], [D_ds; zeros(size(D_ds))], colorResidual, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hResid = plot(s_ds, D_ds, '-', 'Color', colorResidual, 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');

residualMax = prctile(abs(D_residual(edgeMask)), 98) * 1.3;
ylim([-residualMax, residualMax]);
ylabel('D_{2-10} residual (m)', 'FontSize', labelFS, 'Color', colorResidual);
set(gca, 'YColor', colorResidual);

% --- Step 2: LEFT AXIS — white halo + blue depth envelope ---
yyaxis left;

hHalo = plot(s_km(1:ds:end), bands.lowpass.LP_10(1:ds:end), '-', ...
    'Color', [1 1 1], 'LineWidth', 3, 'HandleVisibility', 'off');
hLP = plot(s_km(1:ds:end), bands.lowpass.LP_10(1:ds:end), '-', ...
    'Color', colorLP, 'LineWidth', 2);

ylabel('Depth (m)', 'FontSize', labelFS, 'Color', colorLP);
set(gca, 'YDir', 'reverse', 'YColor', colorLP);

% --- Step 3: Z-order — push green residual behind blue envelope ---
% yyaxis doesn't guarantee left-axis objects render above right-axis.
% Move halo + blue to front of Children, so green sits behind them.
ch = ax1.Children;
frontHandles = [hLP; hHalo];
isFront = ismember(ch, frontHandles);
ax1.Children = [ch(isFront); ch(~isFront)];

% --- Step 4: Samples + labels (drawn AFTER z-reorder, so on top) ---
hSampA = [];
if hasSamples
    sample_LP10 = interp1(s_km, bands.lowpass.LP_10, sample_s_km, 'linear', NaN);
    hSampA = scatter(sample_s_km, sample_LP10, opts.sampleMarkerSize, colorSample, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if opts.showSampleNames
        depthRange = max(bands.lowpass.LP_10(edgeMask)) - min(bands.lowpass.LP_10(edgeMask));
        plotSampleLabels(gca, sample_s_km, sample_LP10, sample_names, depthRange*0.04);
    end
end

% --- Step 5: Furniture ---
drawSegmentLines(segmentBounds, gca, colorBoundary);

text(-0.03, 1.08, 'a', 'Units', 'normalized', 'FontSize', titleFS+4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('Depth profile', 'FontSize', titleFS, 'FontWeight', 'normal');
xlim(xLim);
xlabel('Distance along trench (km)', 'FontSize', labelFS);
set(gca, 'FontSize', tickFS);
grid on;

if hasSamples && ~isempty(hSampA)
    legend([hLP, hResid, hSampA], ...
        {'LP_{>10} depth envelope', 'D_{2-10} residual relief', 'Samples'}, ...
        'Location', 'northeast', 'FontSize', annotFS);
else
    legend([hLP, hResid], ...
        {'LP_{>10} depth envelope', 'D_{2-10} residual relief'}, ...
        'Location', 'northeast', 'FontSize', annotFS);
end
end  % do_a

%% ========================================================================
%  PROFILE B: ROUTING GRADIENTS (dual y-axis)
%  ========================================================================
if do_b
ax2 = nexttile;
hold on;

% --- Left axis: Envelope gradient (LP>10 slope) ---
yyaxis left;
h1 = plot(s_km(1:ds:end), LP10_slope(1:ds:end), '-', ...
    'Color', colorSlopeBroad, 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
slopeMax = prctile(abs(LP10_slope(edgeMask)), 98) * 1.5;
ylim([-slopeMax, slopeMax]);
ylabel('dLP_{>10}/ds (m/km)', 'FontSize', labelFS);
set(gca, 'YColor', colorSlopeBroad);

% --- Right axis: Local variability (median |slope| of D_{2-10}) ---
yyaxis right;
slope_mag_full = interp1(metrics.s_km, metrics.slope_mag, s_km, 'pchip', NaN);
h2 = plot(s_km(1:ds:end), slope_mag_full(1:ds:end), '-', ...
    'Color', colorSlopeLocal, 'LineWidth', 2);
localMax = max(slope_mag_full(~isnan(slope_mag_full))) * 1.1;
ylim([0, localMax]);
ylabel('median|dD_{2-10}/ds| (m/km)', 'FontSize', labelFS);
set(gca, 'YColor', colorSlopeLocal);

% Samples on right axis (local variability) + labels
if hasSamples
    sample_slope = interp1(s_km, slope_mag_full, sample_s_km, 'linear', NaN);
    scatter(sample_s_km, sample_slope, opts.sampleMarkerSize, colorSample, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if opts.showSampleNames
        plotSampleLabels(gca, sample_s_km, sample_slope, sample_names, localMax*0.06, ...
            'OneSided', true);
    end
end

% Segment boundaries
drawSegmentLines(segmentBounds, gca, colorBoundary);

text(-0.03, 1.08, 'b', 'Units', 'normalized', 'FontSize', titleFS+4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('Routing gradients', 'FontSize', titleFS, 'FontWeight', 'normal');
xlim(xLim);
xlabel('Distance along trench (km)', 'FontSize', labelFS);
set(gca, 'FontSize', tickFS);
grid on;

legend([h1, h2], {'dLP_{>10}/ds (envelope gradient)', 'median|dD_{2-10}/ds| (local variability)'}, ...
    'Location', 'northeast', 'FontSize', annotFS);
end  % do_b

%% ========================================================================
%  PROFILE C: SINK SEGMENTATION (Curvature Domains)
%  ========================================================================
if do_c
ax3 = nexttile;
hold on;

% Get curvature data
curv_5_10 = bands.d2D_5_10_mpkm2;

% Smooth for display
smoothWin = max(3, round(5 / ds_km));
curv_smooth = smoothdata(curv_5_10, 'gaussian', smoothWin);

% Y-limits (compute before domain shading)
curvMax = prctile(abs(curv_5_10(edgeMask)), 95) * 1.3;

% --- Identify robust domains and boundaries ---
% Use feature-level stability if available (tests whether the domain as a
% spatial feature persists across parameter perturbations), otherwise fall
% back to boundary-position stability.
robustThreshold = 0.5;
nDom   = length(domains);
nBound = numel(segmentBounds);

hasFeatureStab = isfield(segDiag, 'featureStability') && ...
    ~isempty(segDiag.featureStability);
hasStability = ~isempty(segDiag.boundaryStability);

if hasFeatureStab
    % Feature-level: domain is robust if its feature survives >= threshold
    isDomainRobust = segDiag.featureStability(:) >= robustThreshold;
    % Boundary is robust if both adjacent features are robust
    isRobustBound = true(nBound, 1);
    for b = 1:nBound
        if b <= nDom,     isRobustBound(b) = isRobustBound(b) && isDomainRobust(b); end
        if b+1 <= nDom,   isRobustBound(b) = isRobustBound(b) && isDomainRobust(b+1); end
    end
elseif hasStability
    isRobustBound = segDiag.boundaryStability(:) >= robustThreshold;
    isDomainRobust = true(nDom, 1);
    for d = 1:nDom
        if d > 1,       isDomainRobust(d) = isDomainRobust(d) && isRobustBound(d-1); end
        if d <= nBound,  isDomainRobust(d) = isDomainRobust(d) && isRobustBound(d);   end
    end
else
    isRobustBound  = true(nBound, 1);
    isDomainRobust = true(nDom, 1);
end
nRobustDomains = sum(isDomainRobust);

% --- Domain shading ---
domainAlpha = 0.35;

if ~isempty(domains)
    for d = 1:nDom
        dom = domains(d);
        if dom.meanCurv < 0, fc = colorDomainNeg; else, fc = colorDomainPos; end
        patch([dom.s_start, dom.s_end, dom.s_end, dom.s_start], ...
            [-curvMax, -curvMax, curvMax, curvMax], ...
            fc, 'EdgeColor', 'none', 'FaceAlpha', domainAlpha, ...
            'HandleVisibility', 'off');
    end
end

% Plot smoothed curvature
plot(s_km(1:ds:end), curv_smooth(1:ds:end), '-', ...
    'Color', colorCurvature, 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

ylim([-curvMax, curvMax]);

% Segment boundaries — styled by convergence CR when available, else by robustness
if hasConv && ~isempty(convData) && isfield(convData, 'composite') && numel(convData.composite) == nBound
    % CR-scaled rendering: line width and alpha proportional to composite CR
    %   Map composite CR from [1, crCap] to [lwMin, lwMax] and [alphaMin, alphaMax]
    crCap    = 2.0;
    lwMin    = 0.4;   lwMax = 2.0;
    alphaMin = 0.25;  alphaMax = 0.85;
    for b = 1:nBound
        cr = min(convData.composite(b), crCap);
        t  = (cr - 1) / (crCap - 1);           % normalised [0, 1]
        lw = lwMin + t * (lwMax - lwMin);
        al = alphaMin + t * (alphaMax - alphaMin);
        rgb_cr = colorBoundary * al + (1 - al);  % blend towards white
        xline(gca, segmentBounds(b), '--', 'Color', rgb_cr, ...
            'LineWidth', lw, 'HandleVisibility', 'off');
    end
elseif hasStability
    rgb_faint = colorBoundary * 0.6 + 0.4;
    for b = 1:nBound
        if isRobustBound(b)
            xline(gca, segmentBounds(b), '--', 'Color', rgb_faint, ...
                'LineWidth', 1.0, 'HandleVisibility', 'off');
        else
            rgb_weak = colorBoundary * 0.3 + 0.7;
            xline(gca, segmentBounds(b), ':', 'Color', rgb_weak, ...
                'LineWidth', 0.6, 'HandleVisibility', 'off');
        end
    end
else
    drawSegmentLines(segmentBounds, gca, colorBoundary);
end

% Domain annotations
if ~isempty(domains)
    for d = 1:nDom
        dom = domains(d);
        s_mid = (dom.s_start + dom.s_end) / 2;
        if dom.meanCurv < 0
            label = sprintf('D%d\n\\sigma=%.0fm', d, dom.meanRoughness);
            yPos = -curvMax * 0.7;
        else
            label = sprintf('D%d', d);
            yPos = curvMax * 0.7;
        end
        text(s_mid, yPos, label, 'FontSize', annotFS-1, 'Color', [0 0 0], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% Panel label (placed after axes furniture below)

% Samples
if hasSamples
    sample_curv = interp1(s_km, curv_smooth, sample_s_km, 'linear', NaN);
    hSampC = scatter(sample_s_km, sample_curv, opts.sampleMarkerSize, colorSample, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if opts.showSampleNames
        plotSampleLabels(gca, sample_s_km, sample_curv, sample_names, curvMax*0.08);
    end
end

ylabel('\kappa_{5-10} (m/km^2)', 'FontSize', labelFS);
text(-0.03, 1.08, 'c', 'Units', 'normalized', 'FontSize', titleFS+4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('Sink segmentation', 'FontSize', titleFS, 'FontWeight', 'normal');
xlim(xLim);
xlabel('Distance along trench (km)', 'FontSize', labelFS);
set(gca, 'FontSize', tickFS);
box on;
grid on;

% Mirror tick marks on right y-axis (no label, just numbers + ticks)
yyaxis right;
ylim([-curvMax, curvMax]);
set(gca, 'YColor', [0.15 0.15 0.15], 'FontSize', tickFS);
ylabel('');
yyaxis left;  % return to left axis as active

% Legend — 3 entries: curvature line, trap fill, bypass fill
hCurvLine = plot(NaN, NaN, '-', 'Color', colorCurvature, 'LineWidth', 1.5);
hTrapLeg  = patch(NaN, NaN, colorDomainNeg, 'EdgeColor', 'none', 'FaceAlpha', domainAlpha);
hBypassLeg= patch(NaN, NaN, colorDomainPos, 'EdgeColor', 'none', 'FaceAlpha', domainAlpha);
legendEntries = [hCurvLine, hTrapLeg, hBypassLeg];
legendLabels  = {'\kappa_{5-10} curvature', 'Trapping', 'Bypass'};
if hasSamples && exist('hSampC', 'var')
    legendEntries(end+1) = hSampC;
    legendLabels{end+1} = 'Samples';
end
legend(legendEntries, legendLabels, 'Location', 'northeast', 'FontSize', annotFS);
end  % do_c

%% ========================================================================
%  PROFILE D: SINK STRENGTH (Heatmap Strip)
%  ========================================================================
if do_d
ax4 = nexttile;

% Build 3-row heatmap
heatData = [metrics.rms_norm'; metrics.mad_norm'; metrics.beta_norm'];

% Display
imagesc(metrics.s_km, 1:3, heatData);
colormap(gca, crameri_batlow(64));
caxis([0 1]);

% Y-axis labels (β computed from detrended raw depth, not band-passed)
set(gca, 'YTick', 1:3, 'YTickLabel', {'RMS', 'MAD', '\beta_{raw}'}, 'FontSize', tickFS);
xlim(xLim);

% Segment boundaries (white on heatmap)
hold on;
drawSegmentLines(segmentBounds, gca, [1 1 1]);

% Samples (small markers on top) with labels
if hasSamples
    scatter(sample_s_km, 2*ones(size(sample_s_km)), opts.sampleMarkerSize*0.5, ...
        'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if opts.showSampleNames
        plotSampleLabels(gca, sample_s_km, 2*ones(size(sample_s_km)), sample_names, 0.3, ...
            'OneSided', true);
    end
end

text(-0.03, 1.08, 'd', 'Units', 'normalized', 'FontSize', titleFS+4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('Roughness metrics', 'FontSize', titleFS, 'FontWeight', 'normal');
xlabel('Distance along trench (km)', 'FontSize', labelFS);

% Colorbar
cb = colorbar('Location', 'eastoutside');
ylabel(cb, 'Normalized', 'FontSize', annotFS);
end  % do_d

%% ========================================================================
%  PROFILE E: TRAPPING INDEX (Composite)
%  ========================================================================
if do_e
ax5 = nexttile;
drawTIcurve(ax5, bands, metrics, opts, segmentBounds, ...
    hasSamples, sample_s_km, sample_names, ...
    xLim, colorTI, colorSample, colorBoundary, ...
    titleFS, labelFS, tickFS, annotFS, ds);
text(ax5, -0.03, 1.08, 'e', 'Units', 'normalized', ...
    'FontSize', titleFS+4, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title(ax5, 'Trapping Index', 'FontSize', titleFS, 'FontWeight', 'normal');
end  % do_e

displayLen_km = diff(xLim);  %#ok<NASGU>

%% Link x-axes across rendered panels
allAx = [ax1, ax2, ax3, ax4, ax5];
allAx = allAx(~cellfun(@isempty, num2cell(allAx)));
if numel(allAx) > 1
    linkaxes(allAx, 'x');
end

%% --- Edge degradation overlay ---
% The xLim already excludes fully-NaN regions. But Gaussian filters and
% windowed metrics degrade gradually before going NaN — there is a ramp-up
% zone where values are technically valid but physically meaningless (edge
% artifacts: unrealistic shallowing, TI saturation, curvature spikes).
%
% Detection: scan inward from each edge of xLim along TI_full to find
% where TI drops below a saturation threshold. TI saturates at ~1.0 in
% edge-contaminated zones because all three components (low slope, high
% roughness, concave curvature) are driven to extreme values by filter
% ramp-up. The first point where TI < threshold marks the start of
% trustworthy data.
%
% A single mask applied to all panels is correct here because TI is
% the composite of all underlying metrics — if TI is saturated, ALL
% panels are affected.

edgeSatThreshold = 0.95;   % TI above this at the edge = degraded
TI_disp = metrics.TI_full(:)';
s_disp  = s_km(:)';

% Restrict to displayed range
inView = s_disp >= xLim(1) & s_disp <= xLim(2);
s_view  = s_disp(inView);
TI_view = TI_disp(inView);

% Scan from left: find first index where TI drops below threshold
leftDeg = xLim(1);   % default: no degraded zone
for ii = 1:numel(TI_view)
    if ~isnan(TI_view(ii)) && TI_view(ii) < edgeSatThreshold
        leftDeg = s_view(ii);
        break;
    end
end

% Scan from right: find first index (from right) where TI drops below threshold
rightDeg = xLim(2);  % default: no degraded zone
for ii = numel(TI_view):-1:1
    if ~isnan(TI_view(ii)) && TI_view(ii) < edgeSatThreshold
        rightDeg = s_view(ii);
        break;
    end
end

% Draw gray patches on all rendered panels
panelAxes = allAx;
if leftDeg > xLim(1) || rightDeg < xLim(2)
    fprintf('Edge degradation overlay: left [%.0f to %.0f km], right [%.0f to %.0f km]\n', ...
        xLim(1), leftDeg, rightDeg, xLim(2));
    for iax = 1:numel(panelAxes)
        ax = panelAxes(iax);
        yl = get(ax, 'YLim');
        yP = [yl(1) yl(1) yl(2) yl(2)];
        hold(ax, 'on');
        if leftDeg > xLim(1)
            patch(ax, [xLim(1) leftDeg leftDeg xLim(1)], yP, ...
                [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        if rightDeg < xLim(2)
            patch(ax, [rightDeg xLim(2) xLim(2) rightDeg], yP, ...
                [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        hold(ax, 'off');
    end
end

%% Build results structure
results = struct();
results.metrics = metrics;
results.segmentBounds = segmentBounds;
results.domains = domains;
results.nDomains = length(domains);
if hasConv && ~isempty(convData)
    results.convergence = convData;
end

%% Save if requested
if ~isempty(opts.outputPath)
    [outDir, outName, ~] = fileparts(opts.outputPath);
    if isempty(outDir), outDir = '.'; end
    
    % PNG
    print(fig, fullfile(outDir, [outName, '_1DSink_Master.png']), '-dpng', '-r600');
    
    % PDF
    exportgraphics(fig, fullfile(outDir, [outName, '_1DSink_Master.pdf']), 'ContentType', 'vector');
    
    fprintf('Saved: %s_1DSink_Master.png and .pdf\n', fullfile(outDir, outName));
end

end

%% ========================================================================
% HELPER FUNCTIONS
%% ========================================================================

function [hasSamples, sample_s_km, sample_names, nSamp] = processSamples(samples)
%PROCESSSAMPLES Extract sample positions and names from table/struct
%   Priority for names: OriginalID (from A04) > ID (from A00) > name > SampleID
    hasSamples = ~isempty(samples);
    sample_s_km = [];
    sample_names = {};
    nSamp = 0;
    
    if ~hasSamples
        return;
    end
    
    if istable(samples)
        % Get along-axis position
        if ismember('s_axis', samples.Properties.VariableNames)
            sample_s_km = samples.s_axis / 1000;
        elseif ismember('s_km', samples.Properties.VariableNames)
            sample_s_km = samples.s_km;
        end
        % Get sample names - check OriginalID first (A04 output)
        if ismember('OriginalID', samples.Properties.VariableNames)
            sample_names = samples.OriginalID;
        elseif ismember('ID', samples.Properties.VariableNames)
            sample_names = samples.ID;
        elseif ismember('name', samples.Properties.VariableNames)
            sample_names = samples.name;
        elseif ismember('SampleID', samples.Properties.VariableNames)
            sample_names = samples.SampleID;
        end
    elseif isstruct(samples)
        if isfield(samples, 's_axis')
            sample_s_km = samples.s_axis / 1000;
        elseif isfield(samples, 's_km')
            sample_s_km = samples.s_km;
        end
        if isfield(samples, 'OriginalID')
            sample_names = samples.OriginalID;
        elseif isfield(samples, 'ID')
            sample_names = samples.ID;
        elseif isfield(samples, 'name')
            sample_names = samples.name;
        elseif isfield(samples, 'SampleID')
            sample_names = samples.SampleID;
        end
    end
    
    % Convert to cell array, preserve original names
    sample_names = convertToNameCell(sample_names, nSamp);
    nSamp = numel(sample_s_km);
    hasSamples = nSamp > 0;
end

function nameCell = convertToNameCell(names, ~)
%CONVERTTONAMECELL Convert various name formats to cell array of strings
%   Preserves original names - does NOT convert numeric to "S#" format
    if isempty(names)
        nameCell = {};
    elseif isstring(names)
        nameCell = cellstr(names);
    elseif iscell(names)
        nameCell = names;
    elseif ischar(names)
        nameCell = {names};
    elseif isnumeric(names)
        % Convert numeric to string but preserve original number
        nameCell = arrayfun(@(x) num2str(x), names, 'UniformOutput', false);
    else
        nameCell = {};
    end
end

function plotSampleLabels(ax, s_km, y_vals, names, offset, varargin)
%PLOTSAMPLELABELS Label samples with white-box text and smart collision avoidance
%
%   Uses MATLAB's built-in BackgroundColor on text objects for a clean,
%   perfectly-centered white box behind each label.  Detects clusters of
%   close samples and staggers labels vertically in tiers so they don't
%   overlap.  Optional 'OneSided' mode keeps all labels above the data
%   point (useful when the y-axis starts at 0 or when YDir is reversed).
%
%   Syntax:
%     plotSampleLabels(ax, s_km, y_vals, names, offset)
%     plotSampleLabels(..., 'OneSided', true)
%     plotSampleLabels(..., 'FontSize', 8)
%     plotSampleLabels(..., 'MinSep', 100)
%     plotSampleLabels(..., 'FontColor', [0 0 0])
%     plotSampleLabels(..., 'BgColor',   [1 1 1])
%
%   Inputs:
%     ax       - axes handle
%     s_km     - along-axis positions of samples (km)
%     y_vals   - y-coordinates on current plot where sample dots sit
%     names    - cell array of label strings
%     offset   - base vertical shift from data point (in data units)
%
%   Name-Value Options:
%     'OneSided'  (false) - true = labels only above data points on screen
%     'FontSize'  (7)     - label font size in points
%     'MinSep'    (80)    - samples closer than this (km) get staggered
%     'FontColor' ([0,0,0]) - text color
%     'BgColor'   ([1,1,1]) - background box color

    % ---- Parse optional name-value pairs ----
    p = inputParser;
    addParameter(p, 'OneSided',  false);
    addParameter(p, 'FontSize',  7);
    addParameter(p, 'MinSep',    80);
    addParameter(p, 'FontColor', [0 0 0]);
    addParameter(p, 'BgColor',   [1 1 1]);
    parse(p, varargin{:});
    o = p.Results;

    s_km  = s_km(:);
    y_vals = y_vals(:);
    nSamp = numel(s_km);
    if nSamp == 0, return; end

    isReversed = strcmp(get(ax, 'YDir'), 'reverse');

    % Automatically use OneSided for reversed axes (depth plots) since
    % "below on screen" = deeper and labels would fall off the axis.
    if isReversed
        o.OneSided = true;
    end

    % ---- Tier multipliers for vertical stagger ----
    % Gaps between consecutive tiers are ~1.3× the base offset, ensuring
    % that at 7pt font, adjacent-tier labels don't overlap vertically.
    maxTiers = 6;
    if o.OneSided
        % All tiers go upward on screen (positive for normal, negative for reversed)
        tierMult = [1, 2.3, 3.6, 4.9, 6.2, 7.5];
    else
        % Alternating above / below with wide spacing
        tierMult = [1, -1.3, 2.6, -2.6, 4.0, -4.0];
    end

    % ---- Assign vertical tiers from consecutive proximity ----
    [~, sortIdx] = sort(s_km);
    tier = zeros(nSamp, 1);                % tier in sorted order

    for i = 2:nSamp
        prevPos = s_km(sortIdx(i-1));
        currPos = s_km(sortIdx(i));
        if abs(currPos - prevPos) < o.MinSep
            tier(i) = mod(tier(i-1) + 1, maxTiers);
        else
            tier(i) = 0;                  % reset when gap is wide enough
        end
    end

    % Map tiers back to original sample order
    tierOrig = zeros(nSamp, 1);
    tierOrig(sortIdx) = tier;

    % ---- Get axis limits for clamping ----
    yl = ylim(ax);
    yRange = abs(yl(2) - yl(1));

    % ---- Draw labels + thin leader lines ----
    for i = 1:nSamp
        if isnan(y_vals(i)) || i > numel(names) || isempty(names{i})
            continue;
        end

        % Staggered y-position
        t = tierOrig(i);
        mult = tierMult(mod(t, maxTiers) + 1);

        % For reversed axes, "up on screen" = smaller y-value
        if isReversed
            y_label = y_vals(i) - mult * abs(offset);
        else
            y_label = y_vals(i) + mult * abs(offset);
        end

        % Clamp to axis limits (keep 4% margin from edges)
        margin = yRange * 0.04;
        y_label = max(y_label, yl(1) + margin);
        y_label = min(y_label, yl(2) - margin);

        % Vertical alignment based on position relative to data point
        if (y_label >= y_vals(i) && ~isReversed) || ...
           (y_label <= y_vals(i) &&  isReversed)
            vAlign = 'bottom';     % label sits above on screen
        else
            vAlign = 'top';        % label sits below on screen
        end

        % Draw thin leader line from data point to label (if far enough)
        labelDist = abs(y_label - y_vals(i));
        if labelDist > abs(offset) * 1.5
            plot(ax, [s_km(i) s_km(i)], [y_vals(i) y_label], '-', ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 0.4, ...
                'HandleVisibility', 'off');
        end

        % Single text object with white background box
        text(ax, s_km(i), y_label, names{i}, ...
            'FontSize',            o.FontSize, ...
            'FontWeight',          'bold', ...
            'Color',               o.FontColor, ...
            'BackgroundColor',     o.BgColor, ...
            'Margin',              0.5, ...
            'EdgeColor',           'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   vAlign, ...
            'Clipping',            'on', ...
            'HandleVisibility',    'off');
    end
end

function drawSegmentLines(bounds, ax, rgb)
%DRAWSEGMENTLINES Draw vertical lines at segment boundaries using xline
    if isempty(bounds)
        return;
    end
    
    % Make color slightly faint
    rgb_faint = rgb * 0.6 + 0.4;
    
    for b = bounds(:)'
        xline(ax, b, '--', 'Color', rgb_faint, 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
end

function cmap = crameri_batlow(n)
%CRAMERI_BATLOW  Perceptually uniform, colorblind-safe sequential colormap
%   Attempt at approximate reproduction of Fabio Crameri's "batlow"
%   colormap from key control-point RGB values.
%
%   Properties (following Crameri et al., 2020, Nature Communications):
%     - Perceptually uniform: equal data steps → equal perceived steps
%     - Monotonic lightness:  reads correctly printed in grayscale
%     - Colorblind-safe:      distinguishable under all CVD types
%
%   Usage:
%     colormap(gca, crameri_batlow(64));
%
%   Reference:
%     Crameri, F., Shephard, G.E. & Heron, P.J. (2020). The misuse of
%     colour in science communication. Nature Communications, 11, 5444.
%     https://doi.org/10.1038/s41467-020-19160-7
%
%   For the authoritative version, install the Scientific Colour Maps
%   package from https://www.fabiocrameri.ch/colourmaps/

    if nargin < 1, n = 256; end

    % 17 control points sampled from the published batlow 256-node table.
    % Position (0–1) and RGB triplet.
    nodes = [ ...
        0.000,  0.005,  0.098,  0.350;   % dark indigo
        0.063,  0.024,  0.141,  0.390;
        0.125,  0.055,  0.189,  0.420;
        0.188,  0.098,  0.245,  0.432;
        0.250,  0.153,  0.296,  0.427;   % steel blue
        0.313,  0.212,  0.346,  0.396;
        0.375,  0.280,  0.396,  0.348;   % teal-olive transition
        0.438,  0.358,  0.441,  0.298;
        0.500,  0.448,  0.476,  0.264;   % olive
        0.563,  0.544,  0.507,  0.247;
        0.625,  0.641,  0.530,  0.250;   % ochre
        0.688,  0.736,  0.545,  0.268;
        0.750,  0.820,  0.556,  0.296;   % warm tan
        0.813,  0.886,  0.574,  0.350;
        0.875,  0.935,  0.637,  0.438;   % salmon-pink
        0.938,  0.968,  0.759,  0.522;
        1.000,  0.992,  0.906,  0.144];  % bright yellow

    t    = nodes(:,1);
    rgb  = nodes(:,2:4);

    % Interpolate each channel to n levels
    ti = linspace(0, 1, n)';
    cmap = [ ...
        interp1(t, rgb(:,1), ti, 'pchip'), ...
        interp1(t, rgb(:,2), ti, 'pchip'), ...
        interp1(t, rgb(:,3), ti, 'pchip')];

    % Clamp to valid range
    cmap = max(0, min(1, cmap));
end

