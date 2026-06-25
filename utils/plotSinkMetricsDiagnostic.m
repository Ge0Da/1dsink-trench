function [fig, metrics] = plotSinkMetricsDiagnostic(bands, opts)
%PLOTSINKMETRICSDIAGNOSTIC Diagnostic plot for roughness metrics computation
%
%   fig = plotSinkMetricsDiagnostic(bands)
%   [fig, metrics] = plotSinkMetricsDiagnostic(bands, 'Name', value, ...)
%
%   PURPOSE: Quality control and understanding of sink strength metrics.
%            Shows raw values, distributions, and inter-metric correlations.
%
%   DISPLAY (6 panels):
%     Panel 1: RMS roughness along axis (raw values in meters)
%     Panel 2: MAD roughness along axis (raw values in meters)
%     Panel 3: Spectral slope β along axis (dimensionless)
%     Panel 4: Histograms of all three metrics
%     Panel 5: RMS vs MAD scatter (correlation check)
%     Panel 6: β vs RMS scatter (spectrum-amplitude relationship)
%
%   INPUTS:
%       bands   - Output from computeTrenchBands()
%
%   NAME-VALUE OPTIONS:
%       'title'             - Figure title
%       'windowSize_km'     - Moving window for metrics (default: 50)
%       'samples'           - Sample table/struct (optional)
%       'outputPath'        - Save figure if provided
%
%   OUTPUTS:
%       fig     - Figure handle
%       metrics - Computed sink metrics structure
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

arguments
    bands struct
    opts.title (1,:) char = 'Sink Metrics Diagnostic — RMS, MAD, β Computation'
    opts.windowSize_km (1,1) double = 50
    opts.windowOverlap (1,1) double = 0.8
    opts.w_slope (1,1) double = 1.0
    opts.w_rough (1,1) double = 1.2
    opts.w_curv (1,1) double = 0.8
    opts.normMode (1,:) char {mustBeMember(opts.normMode,{'absolute','percentile'})} = 'absolute'
    opts.slopeLo_mpkm (1,1) double = 0.5
    opts.slopeHi_mpkm (1,1) double = 15.0
    opts.roughLo_m (1,1) double = 2.0
    opts.roughHi_m (1,1) double = 80.0
    opts.curvLo_mpkm2 (1,1) double = -0.01
    opts.curvHi_mpkm2 (1,1) double = 0.05
    opts.pctileLow (1,1) double = 5
    opts.pctileHigh (1,1) double = 95
    opts.betaDetrend (1,:) char {mustBeMember(opts.betaDetrend,{'lp50','lp1000'})} = 'lp50'
    opts.betaFitBand (1,2) double = [0.1, 0.5]
    opts.samples = []
    opts.showSampleNames (1,1) logical = true  % Accepted for interface consistency
    opts.sampleMarkerSize (1,1) double = 40
    opts.outputPath (1,:) char = ''
end

%% Compute metrics
metrics = computeSinkMetrics(bands, ...
    'windowSize_km', opts.windowSize_km, ...
    'windowOverlap', opts.windowOverlap, ...
    'w_slope', opts.w_slope, ...
    'w_rough', opts.w_rough, ...
    'w_curv', opts.w_curv, ...
    'normMode', opts.normMode, ...
    'slopeLo_mpkm', opts.slopeLo_mpkm, ...
    'slopeHi_mpkm', opts.slopeHi_mpkm, ...
    'roughLo_m', opts.roughLo_m, ...
    'roughHi_m', opts.roughHi_m, ...
    'curvLo_mpkm2', opts.curvLo_mpkm2, ...
    'curvHi_mpkm2', opts.curvHi_mpkm2, ...
    'pctileLow', opts.pctileLow, ...
    'pctileHigh', opts.pctileHigh, ...
    'betaDetrend', opts.betaDetrend, ...
    'betaFitBand', opts.betaFitBand);

%% Extract data
s_km = bands.s_km;
xLim = [0, bands.length_km];

%% Process samples
[hasSamples, sample_s_km] = processSamples(opts.samples);

%% Colors
colorRMS = [0.2 0.5 0.8];   % Blue
colorMAD = [0.8 0.4 0.2];   % Orange
colorBeta = [0.5 0.3 0.7];  % Purple
colorSample = [0.95 0.6 0.1];

%% Create figure
fig = figure('Name', 'Sink Metrics Diagnostic', 'Position', [50 50 1500 900], 'Color', 'w');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel 1: RMS along axis
nexttile;
hold on;
plot(metrics.s_km, metrics.rms, '-', 'Color', colorRMS, 'LineWidth', 1.5);
if hasSamples
    sample_rms = interp1(metrics.s_km, metrics.rms, sample_s_km, 'linear', NaN);
    scatter(sample_s_km, sample_rms, 40, colorSample, 'filled', 'MarkerEdgeColor', 'k');
end
ylabel('RMS (m)', 'FontSize', 10);
title('RMS roughness: amplitude of D_{2-10} in windows', 'FontSize', 11, 'FontWeight', 'normal');
text(-0.08, 1.06, 'a', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlim(xLim);
grid on;

% Stats annotation
rms_valid = metrics.rms(~isnan(metrics.rms));
text(xLim(2)*0.02, max(rms_valid)*0.95, ...
    sprintf('Mean: %.1f m\nStd: %.1f m\nRange: %.1f-%.1f m', ...
    mean(rms_valid), std(rms_valid), min(rms_valid), max(rms_valid)), ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'BackgroundColor', 'w');

%% Panel 2: MAD along axis
nexttile;
hold on;
plot(metrics.s_km, metrics.mad, '-', 'Color', colorMAD, 'LineWidth', 1.5);
if hasSamples
    sample_mad = interp1(metrics.s_km, metrics.mad, sample_s_km, 'linear', NaN);
    scatter(sample_s_km, sample_mad, 40, colorSample, 'filled', 'MarkerEdgeColor', 'k');
end
ylabel('MAD (m)', 'FontSize', 10);
title('MAD roughness: robust measure of D_{2-10} variability', 'FontSize', 11, 'FontWeight', 'normal');
text(-0.08, 1.06, 'b', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlim(xLim);
grid on;

% Stats annotation
mad_valid = metrics.mad(~isnan(metrics.mad));
text(xLim(2)*0.02, max(mad_valid)*0.95, ...
    sprintf('Mean: %.1f m\nStd: %.1f m\nRange: %.1f-%.1f m', ...
    mean(mad_valid), std(mad_valid), min(mad_valid), max(mad_valid)), ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'BackgroundColor', 'w');

%% Panel 3: Beta along axis
nexttile;
hold on;
plot(metrics.s_km, metrics.beta, '-', 'Color', colorBeta, 'LineWidth', 1.5);
yline(2, 'k--', 'LineWidth', 0.8);  % Reference: β=2 is Brownian
if hasSamples
    sample_beta = interp1(metrics.s_km, metrics.beta, sample_s_km, 'linear', NaN);
    scatter(sample_s_km, sample_beta, 40, colorSample, 'filled', 'MarkerEdgeColor', 'k');
end
ylabel('\beta (spectral slope)', 'FontSize', 10);
xlabel('Distance along trench axis (km)', 'FontSize', 10);
title('Spectral slope \beta: PSD fit in 2-10 km band (\beta > 2 = smooth, \beta < 2 = rough)', ...
    'FontSize', 11, 'FontWeight', 'normal');
text(-0.08, 1.06, 'c', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlim(xLim);
grid on;

% Stats annotation
beta_valid = metrics.beta(~isnan(metrics.beta));
text(xLim(2)*0.02, max(beta_valid)*0.95, ...
    sprintf('Mean: %.2f\nStd: %.2f\nRange: %.2f-%.2f', ...
    mean(beta_valid), std(beta_valid), min(beta_valid), max(beta_valid)), ...
    'FontSize', 9, 'VerticalAlignment', 'top', 'BackgroundColor', 'w');

%% Panel 4: Histograms
nexttile;
hold on;

% Normalize for overlay
edges_rms = linspace(min(rms_valid), max(rms_valid), 30);
edges_mad = linspace(min(mad_valid), max(mad_valid), 30);
edges_beta = linspace(min(beta_valid), max(beta_valid), 30);

h1 = histogram(rms_valid, edges_rms, 'Normalization', 'probability', ...
    'FaceColor', colorRMS, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = histogram(mad_valid, edges_mad, 'Normalization', 'probability', ...
    'FaceColor', colorMAD, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

ylabel('Probability', 'FontSize', 10);
xlabel('Value (m)', 'FontSize', 10);
title('Distribution of RMS and MAD', 'FontSize', 11, 'FontWeight', 'normal');
text(-0.10, 1.06, 'd', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend([h1, h2], {'RMS', 'MAD'}, 'Location', 'northeast');
grid on;

% RMS/MAD ratio annotation
ratio = mean(rms_valid) / mean(mad_valid);
text(0.98, 0.95, sprintf('RMS/MAD ratio: %.2f', ratio), ...
    'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'BackgroundColor', 'w');

%% Panel 5: RMS vs MAD scatter
nexttile;
hold on;

scatter(metrics.rms, metrics.mad, 30, colorRMS, 'filled', 'MarkerFaceAlpha', 0.5);

% Fit line
validIdx = ~isnan(metrics.rms) & ~isnan(metrics.mad);
if sum(validIdx) > 5
    p = polyfit(metrics.rms(validIdx), metrics.mad(validIdx), 1);
    xFit = linspace(min(metrics.rms(validIdx)), max(metrics.rms(validIdx)), 100);
    yFit = polyval(p, xFit);
    plot(xFit, yFit, 'k-', 'LineWidth', 1.5);
    
    % Correlation
    r = corr(metrics.rms(validIdx), metrics.mad(validIdx));
    text(0.05, 0.95, sprintf('r = %.3f', r), 'Units', 'normalized', ...
        'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end

% 1:1 line for reference
maxVal = max([max(metrics.rms), max(metrics.mad)]);
plot([0 maxVal], [0 maxVal], 'k:', 'LineWidth', 0.8);

xlabel('RMS (m)', 'FontSize', 10);
ylabel('MAD (m)', 'FontSize', 10);
title('RMS vs MAD: correlation check', 'FontSize', 11, 'FontWeight', 'normal');
text(-0.10, 1.06, 'e', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
axis equal;
grid on;

% Samples
if hasSamples
    scatter(sample_rms, sample_mad, 60, colorSample, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

%% Panel 6: Beta vs RMS scatter
nexttile;
hold on;

scatter(metrics.rms, metrics.beta, 30, colorBeta, 'filled', 'MarkerFaceAlpha', 0.5);

% Reference line at β=2
yline(2, 'k--', 'LineWidth', 0.8);

% Correlation
validIdx = ~isnan(metrics.rms) & ~isnan(metrics.beta);
if sum(validIdx) > 5
    r = corr(metrics.rms(validIdx), metrics.beta(validIdx));
    text(0.05, 0.95, sprintf('r = %.3f', r), 'Units', 'normalized', ...
        'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end

xlabel('RMS (m)', 'FontSize', 10);
ylabel('\beta (spectral slope)', 'FontSize', 10);
title('\beta vs RMS: spectrum-amplitude relationship', 'FontSize', 11, 'FontWeight', 'normal');
text(-0.10, 1.06, 'f', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Samples
if hasSamples
    scatter(sample_rms, sample_beta, 60, colorSample, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

% Interpretation annotation
text(0.98, 0.05, 'High β + High RMS = smooth, large-amplitude relief', ...
    'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'FontAngle', 'italic');
text(0.98, 0.12, 'Low β + High RMS = rough, chaotic relief', ...
    'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'FontAngle', 'italic');

%% Overall title removed for paper publication.
%  Per-panel titles + the a-f panel-letter labels supply identification;
%  the supplementary caption carries the descriptive header.

%% Save if requested
if ~isempty(opts.outputPath)
    [outDir, outName, ~] = fileparts(opts.outputPath);
    if isempty(outDir), outDir = '.'; end
    % Use saveFigure for consistent export with transparency support
    saveFigure(fig, outDir, [outName, '_SinkMetricsDiag'], {'png', 'pdf'}, ...
        'resolution', 300, 'transparent', true, 'verbose', true);
end

end

%% ========================================================================
% HELPER FUNCTIONS
%% ========================================================================

function [hasSamples, sample_s_km] = processSamples(samples)
    hasSamples = ~isempty(samples);
    sample_s_km = [];
    if ~hasSamples, return; end
    
    if istable(samples)
        if ismember('s_axis', samples.Properties.VariableNames)
            sample_s_km = samples.s_axis / 1000;
        elseif ismember('s_km', samples.Properties.VariableNames)
            sample_s_km = samples.s_km;
        end
    elseif isstruct(samples)
        if isfield(samples, 's_axis')
            sample_s_km = samples.s_axis / 1000;
        elseif isfield(samples, 's_km')
            sample_s_km = samples.s_km;
        end
    end
    hasSamples = ~isempty(sample_s_km);
end
