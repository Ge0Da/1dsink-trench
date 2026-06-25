% paper_figures.m
% =========================================================================
%  Build the published figure set from cached pipeline outputs.
% =========================================================================
%
%  Produces the five paper PNGs:
%       output/mainFig2_morphometric.png
%       output/mainFig3_TI_and_correlations.png
%       output/mainFig4a_TI_simplex.png
%       output/suppFigS5_diagnostic_metrics.png
%       output/suppFigS6_calibration_diagnostics.png
%
%  Required inputs:
%       output/1DSink_MasterFigure_results.mat   (from run_master_figure)
%       output/TI_weight_calibration.mat         (from calibrate_TI_weights)
%       data/kamchatka_samples.xlsx              (in-repo sample table)
%
%  Usage:
%       cd 101_s2s-along-trench-v1.0
%       paper_figures
%
%  Author: David Fernandez-Blanco
%  Date:   2026-05-20
% =========================================================================

if ~exist('tTotal', 'var'), clear; clc; end
projRoot = fileparts(fileparts(mfilename('fullpath')));
if isempty(projRoot), projRoot = fileparts(pwd); end
addpath(genpath(projRoot));
cd(projRoot);

fprintf('\n');
fprintf('================================================================\n');
fprintf('  paper_figures.m   build the published figure set\n');
fprintf('================================================================\n\n');

mat_master = 'output/1DSink_MasterFigure_results.mat';
mat_calib  = 'output/TI_weight_calibration.mat';
xlsx_samp  = 'data/kamchatka_samples.xlsx';
for f = {mat_master, mat_calib, xlsx_samp}
    if ~isfile(f{1})
        error('paper_figures:missingInput', '%s not found.', f{1});
    end
end

ld    = load(mat_master);
calib = load(mat_calib);
T     = readtable(xlsx_samp, 'VariableNamingRule', 'preserve');

dpi   = 600;
outDir = 'output';
if ~exist(outDir, 'dir'), mkdir(outDir); end

data = compute_correlations(ld, T);

build_figure_2(ld, dpi, outDir);
build_figure_3(ld, data, dpi, outDir);
build_figure_4a(calib.results, dpi, outDir);
build_figure_S5(ld, dpi, outDir);
build_figure_S6(calib.results, dpi, outDir);

fprintf('\n================================================================\n');
fprintf('  Complete.\n');
fprintf('================================================================\n\n');


%% ========================================================================
%  Data
%  ========================================================================
function data = compute_correlations(ld, T)
%COMPUTE_CORRELATIONS Spearman correlations between TI and pollutant classes.

    mr            = ld.masterResults;
    bands_s_km    = ld.bands.s_km(:);
    TI_full       = mr.metrics.TI_full(:);
    sample_s_km   = ld.projectedSamples.s_km(:);
    nSamples      = numel(sample_s_km);

    validTI    = ~isnan(TI_full) & ~isnan(bands_s_km);
    [s_sorted, order] = sort(bands_s_km(validTI));
    TI_sorted  = TI_full(validTI);
    TI_sorted  = TI_sorted(order);
    TI_interp  = interp1(s_sorted, TI_sorted, sample_s_km, 'linear', NaN);

    colNames = T.Properties.VariableNames;
    Total_OCP = sumOrNaN(table2array(T(:, startsWith(colNames, 'OCP_'))));
    Total_PAH = sumOrNaN(table2array(T(:, startsWith(colNames, 'PAH_'))));

    hydAll       = startsWith(colNames, 'HYD_');
    isD50        = contains(colNames, 'D50_IS');
    isBranched   = contains(colNames, '3_methyl');
    isIsoprenoid = contains(colNames, 'Pristan') | contains(colNames, 'Phytane');
    isOlefin     = contains(colNames, 'hexadecene') | contains(colNames, 'octadecene');
    strictCols   = hydAll & ~isD50 & ~isBranched & ~isIsoprenoid & ~isOlefin;
    Total_HYD    = sumOrNaN(table2array(T(:, strictCols)));

    HM_norm_mean = targetMean(T, colNames);   % calibration target (calib_config)

    [rho_HM, p_HM]   = spearman(TI_interp, HM_norm_mean);
    [rho_PAH, p_PAH] = spearman(TI_interp, Total_PAH);
    [rho_OCP, p_OCP] = spearman(TI_interp, Total_OCP);
    [rho_HYD, p_HYD] = spearman(TI_interp, Total_HYD);

    data = struct();
    data.TI_interp     = TI_interp;
    data.sample_s_km   = sample_s_km;
    data.nSamples      = nSamples;
    data.HM_norm_mean  = HM_norm_mean;
    data.Total_PAH     = Total_PAH;
    data.Total_OCP     = Total_OCP;
    data.Total_HYD     = Total_HYD;
    data.rho_HM        = rho_HM;   data.p_HM  = p_HM;
    data.rho_PAH       = rho_PAH;  data.p_PAH = p_PAH;
    data.rho_OCP       = rho_OCP;  data.p_OCP = p_OCP;
    data.rho_HYD       = rho_HYD;  data.p_HYD = p_HYD;

    function v = getCol(Tbl, colname)
        if ~any(strcmp(Tbl.Properties.VariableNames, colname))
            v = nan(height(Tbl), 1);
            return
        end
        raw = Tbl.(colname);
        if iscell(raw)
            v = nan(numel(raw), 1);
            for jj = 1:numel(raw)
                cell_j = raw{jj};
                if isnumeric(cell_j) && isscalar(cell_j)
                    v(jj) = double(cell_j);
                elseif ischar(cell_j) || isstring(cell_j)
                    v(jj) = str2double(cell_j);
                end
            end
        else
            v = double(raw);
        end
    end
end

function v = sumOrNaN(M)
%SUMORNAN Column-wise omit-NaN sum that returns NaN for fully-NaN rows.
    if isempty(M)
        v = nan(size(M, 1), 1); return
    end
    v = sum(M, 2, 'omitnan');
    v(all(isnan(M), 2)) = NaN;
end

function [rho, p] = spearman(x, y)
    valid = ~isnan(x) & ~isnan(y);
    if sum(valid) < 4
        rho = NaN; p = NaN; return
    end
    [rho, p] = corr(x(valid), y(valid), 'type', 'Spearman');
end


%% ========================================================================
%  Figure 2 -- Morphometric decomposition (profiles a, b, c, d)
%  ========================================================================
function build_figure_2(ld, dpi, outDir)
%BUILD_FIGURE_2 Render Fig. 2 using plotTrenchBands_1DSink with 4 profiles.

    fprintf('Building Fig. 2 (morphometric, profiles a-d) ...\n');
    cfg = ld.cfg;
    fig = plotTrenchBands_1DSink(ld.bands, ...
        'panels',                {'a','b','c','d'}, ...
        'samples',               ld.projectedSamples, ...
        'showSampleNames',       cfg.samples.showNames, ...
        'sampleMarkerSize',      cfg.plot.sampleMarkerSize, ...
        'downsamplePlot',        cfg.plot.downsampleFactor, ...
        'windowSize_km',         cfg.metrics.windowSize_km, ...
        'windowOverlap',         cfg.metrics.windowOverlap, ...
        'w_slope',               cfg.metrics.w_slope, ...
        'w_rough',               cfg.metrics.w_rough, ...
        'w_curv',                cfg.metrics.w_curv, ...
        'normMode',              cfg.metrics.normMode, ...
        'slopeLo_mpkm',          cfg.metrics.slopeLo_mpkm, ...
        'slopeHi_mpkm',          cfg.metrics.slopeHi_mpkm, ...
        'roughLo_m',             cfg.metrics.roughLo_m, ...
        'roughHi_m',             cfg.metrics.roughHi_m, ...
        'curvLo_mpkm2',          cfg.metrics.curvLo_mpkm2, ...
        'curvHi_mpkm2',          cfg.metrics.curvHi_mpkm2, ...
        'pctileLow',             cfg.metrics.pctileLow, ...
        'pctileHigh',            cfg.metrics.pctileHigh, ...
        'betaDetrend',           cfg.metrics.betaDetrend, ...
        'betaFitBand',           cfg.metrics.betaFitBand_cyc_km, ...
        'minDomainLength_km',    cfg.segmentation.minDomainLength_km, ...
        'persistenceLength_km',  cfg.segmentation.persistenceLength_km, ...
        'deadbandMode',          cfg.segmentation.deadbandMode, ...
        'deadbandAmplitude',     cfg.segmentation.deadbandAmplitude_fixed, ...
        'deadbandFactor',        cfg.segmentation.deadbandFactor_percentile, ...
        'deadbandSigmas',        cfg.segmentation.deadbandSigmas_mad, ...
        'deadbandWindow_km',     cfg.segmentation.deadbandWindow_km, ...
        'confidenceWeights',     cfg.segmentation.confidenceWeights, ...
        'dwellScale_km',         cfg.segmentation.dwellScale_km, ...
        'localAmplitudeRadius_km', cfg.segmentation.localAmplitudeRadius_km, ...
        'maxBoundaries',         cfg.segmentation.maxBoundaries, ...
        'segmentSmoothWindow_km', cfg.segmentation.smoothWindow_km, ...
        'stabilityTest',         cfg.segmentation.stabilityTest, ...
        'stabilityPerturbation', cfg.segmentation.stabilityPerturbation, ...
        'stabilityGrid',         cfg.segmentation.stabilityGrid, ...
        'stabilitySnapRadius_km',cfg.segmentation.stabilitySnapRadius_km, ...
        'minStability',          cfg.segmentation.minStability, ...
        'precomputedSegmentation', ld.precomputedSeg);
    save_figure(fig, fullfile(outDir, 'mainFig2_morphometric.png'), dpi);
end


%% ========================================================================
%  Figure 3 -- TI banner (a) above 4 contaminant scatters (b i-iv)
%  ========================================================================
function build_figure_3(ld, data, dpi, outDir)
%BUILD_FIGURE_3 Compose a single figure: TI banner + 4 scatters in a row.

    fprintf('Building Fig. 3 (TI banner + 4 scatters) ...\n');

    mr        = ld.masterResults;
    metrics   = mr.metrics;
    bands     = ld.bands;
    samples   = ld.projectedSamples;

    % Display window from edgeMask + valid TI region.
    s_km      = bands.s_km;
    edgeMask  = bands.edgeMask;
    validDisp = edgeMask(:) & ~isnan(metrics.TI_full(:));
    validIdx  = find(validDisp);
    if isempty(validIdx)
        xLim = [s_km(1), s_km(end)];
    else
        xLim = [s_km(validIdx(1)), s_km(validIdx(end))];
    end

    % Segment boundaries (clip to display range, if precomputed).
    if isfield(ld, 'precomputedSeg') && isfield(ld.precomputedSeg, 'boundaries')
        allBounds     = ld.precomputedSeg.boundaries;
        segmentBounds = allBounds(allBounds >= xLim(1) & allBounds <= xLim(2));
    else
        segmentBounds = [];
    end

    fig = figure('Units', 'centimeters', 'Position', [2 2 40 13], ...
                 'Color', 'w', 'Visible', 'off');

    % Layout: banner across the top; 4 scatter panels in a row below.
    fig_left   = 0.05;
    fig_right  = 0.95;
    panel_gap  = 0.030;            % tightened from 0.045
    n_scatter  = 4;
    panel_w    = (fig_right - fig_left - (n_scatter - 1) * panel_gap) / n_scatter;
    panel_lefts = fig_left + (0:n_scatter-1) * (panel_w + panel_gap);

    % Vertical layout: banner at top (with title above and its own x-label
    % below), scatter row at bottom.
    banner_y   = 0.69;
    banner_h   = 0.20;
    scatter_y  = 0.12;
    scatter_h  = 0.38;

    % Panel-letter labels: both at the same figure-x for visual alignment.
    panel_letter_x = 0.005;
    letter_fs      = 16;

    % --- Panel a: TI banner ---
    axA = axes('Position', [fig_left, banner_y, fig_right - fig_left, banner_h]); %#ok<LAXES>
    drawTIcurve(axA, bands, metrics, ...
        struct('sampleMarkerSize', 50, 'showSampleNames', true), ...
        segmentBounds, ~isempty(samples.s_km), ...
        samples.s_km(:), samples.name, ...
        xLim, [0.9 0.5 0.1], [0.95 0.6 0.1], [0.7 0.2 0.2], ...
        11, 10, 9, 8, 5);
    ht = title(axA, 'Trapping Index (TI)', 'FontSize', 11, 'FontWeight', 'normal');
    set(ht, 'Units', 'normalized', 'Position', [0.5, 1.07, 0], ...
            'VerticalAlignment', 'bottom');
    annotation(fig, 'textbox', ...
        [panel_letter_x, banner_y + banner_h + 0.02, 0.04, 0.05], ...
        'String', 'a', 'EdgeColor', 'none', ...
        'FontSize', letter_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    % --- Panel b: 4 contaminant scatters in a row ---
    panelDef = { ...
        data.HM_norm_mean, 'Mean normalised heavy-metal content',   data.rho_HM,  data.p_HM,  '(i)';
        data.Total_PAH,    'Total polycyclic aromatic hydrocarbons', data.rho_PAH, data.p_PAH, '(ii)';
        data.Total_OCP,    'Total organochlorine pesticides',        data.rho_OCP, data.p_OCP, '(iii)';
        data.Total_HYD,    'Total n-alkanes (C14 to C29)',           data.rho_HYD, data.p_HYD, '(iv)'};

    markerClr = [0.18 0.42 0.78];
    for sp = 1:n_scatter
        ax   = axes('Position', [panel_lefts(sp), scatter_y, panel_w, scatter_h]); %#ok<LAXES>
        yi   = panelDef{sp, 1};
        ylab = panelDef{sp, 2};
        rho  = panelDef{sp, 3};
        pval = panelDef{sp, 4};
        lbl  = panelDef{sp, 5};
        draw_scatter_panel(ax, data.TI_interp, yi, ylab, rho, pval, lbl, markerClr);
    end

    annotation(fig, 'textbox', ...
        [panel_letter_x, scatter_y + scatter_h + 0.02, 0.04, 0.05], ...
        'String', 'b', 'EdgeColor', 'none', ...
        'FontSize', letter_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    save_figure(fig, fullfile(outDir, 'mainFig3_TI_and_correlations.png'), dpi);
end

function draw_scatter_panel(ax, x, y, ylab, rho, pval, label, markerClr)
%DRAW_SCATTER_PANEL One contaminant-vs-TI scatter with linear trend + annotation.
    valid = ~isnan(x) & ~isnan(y);

    scatter(ax, x(valid), y(valid), 45, 'filled', ...
        'MarkerFaceColor', markerClr, 'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeColor', 'none');

    if sum(valid) >= 5
        pf = polyfit(x(valid), y(valid), 1);
        xr = linspace(min(x(valid)), max(x(valid)), 100);
        hold(ax, 'on');
        plot(ax, xr, polyval(pf, xr), 'r-', 'LineWidth', 1.5);
        hold(ax, 'off');
    end

    xlabel(ax, 'Trapping Index (TI)', 'FontSize', 9);
    ylabel(ax, ylab, 'FontSize', 9);
    ht = title(ax, sprintf('%s  %s', label, ylab), ...
        'FontSize', 10, 'FontWeight', 'bold');
    set(ht, 'Units', 'normalized', 'Position', [0.5, 1.07, 0], ...
            'VerticalAlignment', 'bottom');

    if ~isnan(rho)
        if     pval < 0.01, sig = '**';
        elseif pval < 0.05, sig = '*';
        else,               sig = '';
        end
        annStr = sprintf('\\rho = %.3f  p = %.3f%s  n = %d', ...
            rho, pval, sig, sum(valid));
        text(ax, 0.05, 0.95, annStr, ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Color', [0.15 0.15 0.15], ...
            'BackgroundColor', [1 1 1 0.6]);
    end

    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.3);
end


%% ========================================================================
%  Figure 4a -- TI weight simplex (equilateral ternary)
%  ========================================================================
function build_figure_4a(results, dpi, outDir)
%BUILD_FIGURE_4A Ternary contour of Spearman rho on the weight simplex.

    fprintf('Building Fig. 4a (ternary TI weight simplex) ...\n');

    rhoGrid = results.rhoGrid;
    wvals   = results.wvals;
    bestWs  = results.optimal_weights;
    wCur    = results.current_weights_normalised;
    bootW   = results.bootstrap_weights;
    nBoot   = size(bootW, 1);

    V_slope = [0,   0];           % LowSlope vertex (bottom-left)
    V_rough = [1,   0];           % Roughness vertex (bottom-right)
    V_curv  = [0.5, sqrt(3)/2];   % Concavity vertex (top)
    bc2xy   = @(ws, wr, wc) deal( ...
        ws*V_slope(1) + wr*V_rough(1) + wc*V_curv(1), ...
        ws*V_slope(2) + wr*V_rough(2) + wc*V_curv(2));

    [Ws, Wr] = meshgrid(wvals, wvals);
    Ws = Ws';  Wr = Wr';
    Wc = 1 - Ws - Wr;
    valid = Wc >= -1e-10 & ~isnan(rhoGrid);
    [Xv, Yv] = bc2xy(Ws(valid), Wr(valid), max(0, Wc(valid)));
    rho_v = rhoGrid(valid);

    nFine = 240;
    [Xg, Yg] = meshgrid(linspace(-0.02, 1.02, nFine), ...
                        linspace(-0.02, sqrt(3)/2 + 0.02, nFine));
    F = scatteredInterpolant(Xv, Yv, rho_v, 'linear', 'none');
    Rg = F(Xg, Yg);
    Wc_g = 2*Yg / sqrt(3);
    Wr_g = Xg - Yg / sqrt(3);
    Ws_g = 1 - Wr_g - Wc_g;
    Rg(~((Ws_g >= 0) & (Wr_g >= 0) & (Wc_g >= 0))) = NaN;

    fig = figure('Units', 'centimeters', 'Position', [2 2 20 16], ...
                 'Color', 'w', 'Visible', 'off');
    ax = axes('Position', [0.05 0.08 0.65 0.84]); %#ok<LAXES>
    hold(ax, 'on');

    contourf(ax, Xg, Yg, Rg, 20, 'LineColor', 'none');
    colormap(ax, parula(64));
    cb = colorbar(ax, 'Location', 'eastoutside');
    cb.Position = [0.82 0.18 0.035 0.65];
    ylabel(cb, 'Spearman \rho', 'FontSize', 11);

    plot(ax, [V_slope(1) V_rough(1) V_curv(1) V_slope(1)], ...
             [V_slope(2) V_rough(2) V_curv(2) V_slope(2)], 'k-', 'LineWidth', 1.2);

    % Edge tick marks + numeric labels at 20% intervals. Convention: each
    % edge labels the component that INCREASES toward the nearer vertex of
    % that component (standard petrology ternary).
    %   Bottom edge (V_slope -> V_rough):  labels w_rough, 0 -> 1 left to right
    %   Right  edge (V_rough -> V_curv):   labels w_curv,  0 -> 1 bottom to top
    %   Left   edge (V_curv  -> V_slope):  labels w_slope, 0 -> 1 top to bottom
    tickFrac = 0.2:0.2:0.8;
    tickLen  = 0.020;
    tickLabelOffset = 0.012;
    for f = tickFrac
        % Bottom edge -- w_rough = f
        pB        = V_slope + f * (V_rough - V_slope);
        tickB     = pB - tickLen * [0, 1];
        plot(ax, [pB(1) tickB(1)], [pB(2) tickB(2)], 'k-', 'LineWidth', 0.8);
        text(ax, tickB(1), tickB(2) - tickLabelOffset, sprintf('%.1f', f), ...
            'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

        % Right edge -- w_curv = f
        pR        = V_rough + f * (V_curv - V_rough);
        dirR      = [cosd(30), sind(30)];   % outward-perpendicular to right edge
        tickR     = pR + tickLen * dirR;
        plot(ax, [pR(1) tickR(1)], [pR(2) tickR(2)], 'k-', 'LineWidth', 0.8);
        text(ax, tickR(1) + tickLabelOffset*cosd(30), tickR(2) + tickLabelOffset*sind(30), ...
            sprintf('%.1f', f), ...
            'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

        % Left edge -- w_slope = f
        pL        = V_curv + f * (V_slope - V_curv);
        dirL      = [cosd(180-30), sind(180-30)]; % outward-perpendicular to left edge
        tickL     = pL + tickLen * dirL;
        plot(ax, [pL(1) tickL(1)], [pL(2) tickL(2)], 'k-', 'LineWidth', 0.8);
        text(ax, tickL(1) + tickLabelOffset*cosd(180-30), tickL(2) + tickLabelOffset*sind(180-30), ...
            sprintf('%.1f', f), ...
            'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end

    text(ax, V_slope(1)-0.04, V_slope(2)-0.06, 'LowSlope', ...
        'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(ax, V_rough(1)+0.04, V_rough(2)-0.06, 'Roughness', ...
        'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(ax, V_curv(1), V_curv(2)+0.04, 'Concavity', ...
        'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    % 95% bootstrap envelope (convex hull of the central 95% of resampled
    % optimal weights projected into the ternary).
    bX = zeros(nBoot, 1);  bY = zeros(nBoot, 1);
    for k = 1:nBoot
        [bX(k), bY(k)] = bc2xy(bootW(k, 1), bootW(k, 2), bootW(k, 3));
    end
    distSq = (bX - mean(bX)).^2 + (bY - mean(bY)).^2;
    inner  = distSq <= prctile(distSq, 95);
    hEnv = [];
    if sum(inner) >= 3
        iX = bX(inner);  iY = bY(inner);
        K  = convhull(iX, iY);
        hEnv = plot(ax, iX(K), iY(K), 'k--', 'LineWidth', 1.0);
    end

    [xOpt, yOpt] = bc2xy(bestWs(1), bestWs(2), bestWs(3));
    hOpt = plot(ax, xOpt, yOpt, 'rp', 'MarkerSize', 18, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    [xOp, yOp] = bc2xy(wCur(1), wCur(2), wCur(3));
    hAdopt = plot(ax, xOp, yOp, 'ws', 'MarkerSize', 14, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    axis(ax, 'equal');
    axis(ax, 'off');
    set(ax, 'XLim', [-0.20, 1.20], 'YLim', [-0.15, sqrt(3)/2 + 0.10]);

    if ~isempty(hEnv)
        leg = legend(ax, [hOpt, hAdopt, hEnv], ...
            {sprintf('Optimum (\\rho = %.3f)', results.optimal_rho), ...
             sprintf('Adopted weights (\\rho = %.3f)', results.current_rho), ...
             '95% bootstrap envelope'}, ...
            'FontSize', 9, 'Box', 'on');
    else
        leg = legend(ax, [hOpt, hAdopt], ...
            {sprintf('Optimum (\\rho = %.3f)', results.optimal_rho), ...
             sprintf('Adopted weights (\\rho = %.3f)', results.current_rho)}, ...
            'FontSize', 9, 'Box', 'on');
    end
    leg.Color    = [1 1 1];
    % Place the legend in the empty top-right region of the figure, clear of
    % the LowSlope vertex label and the colorbar.
    leg.Units    = 'normalized';
    leg.Position = [0.55 0.78 0.22 0.10];

    hold(ax, 'off');

    save_figure(fig, fullfile(outDir, 'mainFig4a_TI_simplex.png'), dpi);
end


%% ========================================================================
%  Supp Fig. S5 -- RMS / MAD / beta diagnostic
%  ========================================================================
function build_figure_S5(ld, dpi, outDir)
%BUILD_FIGURE_S5 Re-render the 6-panel diagnostic from cached bands.

    fprintf('Building Supp. Fig. S5 (RMS / MAD / beta diagnostic) ...\n');
    cfg = ld.cfg;
    fig = plotSinkMetricsDiagnostic(ld.bands, ...
        'title',            '', ...
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
        'samples',          ld.projectedSamples, ...
        'showSampleNames',  cfg.samples.showNames, ...
        'sampleMarkerSize', cfg.plot.sampleMarkerSize_diagnostic);
    save_figure(fig, fullfile(outDir, 'suppFigS5_diagnostic_metrics.png'), dpi);
end


%% ========================================================================
%  Supp Fig. S6 -- TI calibration diagnostics (a, b, c)
%  ========================================================================
function build_figure_S6(results, dpi, outDir)
%BUILD_FIGURE_S6 Three calibration diagnostic panels: single-component
%  Spearman bar, LOOCV scatter, bootstrap weight distributions.

    fprintf('Building Supp. Fig. S6 (calibration diagnostics, a-c) ...\n');

    fig = figure('Units', 'centimeters', 'Position', [2 2 36 12], ...
                 'Color', 'w', 'Visible', 'off');

    % --- a: single-component Spearman bar chart ---
    axA = axes('Position', [0.06 0.22 0.25 0.66]); %#ok<LAXES>
    compNames = {'LowSlope', 'Roughness', 'Concavity'};
    compRho   = results.component_rho(:);
    compP     = results.component_p(:);
    bar(axA, 1:3, compRho, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
    hold(axA, 'on');
    yline(axA, 0, 'k-', 'LineWidth', 0.5);
    for i = 1:3
        if     compP(i) < 0.01, star = '**';
        elseif compP(i) < 0.05, star = '*';
        else,                   star = 'ns';
        end
        ypos = compRho(i) + sign(compRho(i)) * 0.03;
        text(axA, i, ypos, star, 'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    set(axA, 'XTick', 1:3, 'XTickLabel', compNames, 'FontSize', 9);
    ylabel(axA, 'Spearman \rho vs HM', 'FontSize', 10);
    ylim(axA, [-0.5, 1.0]);
    box(axA, 'on'); grid(axA, 'on');
    title(axA, 'a. Single-component diagnostics', 'FontSize', 11, 'FontWeight', 'normal');

    % --- b: LOOCV scatter ---
    axB = axes('Position', [0.39 0.22 0.25 0.66]); %#ok<LAXES>
    TI_loocv  = results.loocv_TI(:);
    HM_norm   = results.HM_norm_mean(:);
    scatter(axB, TI_loocv, HM_norm, 60, 'filled', ...
        'MarkerFaceColor', [0.8 0.3 0.2], 'MarkerFaceAlpha', 0.8);
    hold(axB, 'on');
    valid = ~isnan(TI_loocv) & ~isnan(HM_norm);
    if sum(valid) >= 2
        pf = polyfit(TI_loocv(valid), HM_norm(valid), 1);
        xr = linspace(min(TI_loocv(valid)), max(TI_loocv(valid)), 100);
        plot(axB, xr, polyval(pf, xr), 'r-', 'LineWidth', 1.5);
    end
    xlabel(axB, 'LOOCV-predicted TI', 'FontSize', 10);
    ylabel(axB, 'Mean HM (normalised)', 'FontSize', 10);
    annStr = sprintf('LOOCV \\rho = %.3f\nFull-data \\rho = %.3f\nShrinkage = %.1f%%', ...
        results.loocv_rho, results.optimal_rho, ...
        100 * results.shrinkage / abs(results.optimal_rho));
    text(axB, 0.05, 0.95, annStr, ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, ...
        'BackgroundColor', [1 1 1 0.7]);
    box(axB, 'on'); grid(axB, 'on');
    title(axB, sprintf('b. LOOCV validation (\\rho = %.3f)', results.loocv_rho), ...
        'FontSize', 11, 'FontWeight', 'normal');

    % --- c: bootstrap weight distributions ---
    axC = axes('Position', [0.72 0.22 0.25 0.66]); %#ok<LAXES>
    bootW   = results.bootstrap_weights;
    medW    = results.bootstrap_weights_median;
    colours = [0.2 0.5 0.8; 0.8 0.3 0.2; 0.3 0.7 0.3];
    hold(axC, 'on');
    for i = 1:3
        [counts, edges] = histcounts(bootW(:, i), 20, 'Normalization', 'probability');
        centres = (edges(1:end-1) + edges(2:end)) / 2;
        bar(axC, centres, counts, 'FaceColor', colours(i, :), 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none', 'BarWidth', 1);
    end
    legend(axC, {'w_{slope}', 'w_{rough}', 'w_{curv}'}, ...
        'Location', 'northeast', 'FontSize', 9);
    for i = 1:3
        xline(axC, medW(i), '--', 'Color', colours(i, :), 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
    xlabel(axC, 'Weight value', 'FontSize', 10);
    ylabel(axC, 'Probability', 'FontSize', 10);
    box(axC, 'on'); grid(axC, 'on');
    nBoot = size(bootW, 1);
    title(axC, sprintf('c. Bootstrap weight distributions (n = %d)', nBoot), ...
        'FontSize', 11, 'FontWeight', 'normal');

    save_figure(fig, fullfile(outDir, 'suppFigS6_calibration_diagnostics.png'), dpi);
end


%% ========================================================================
%  Save helper
%  ========================================================================
function save_figure(fig, path, dpi)
%SAVE_FIGURE Render figure with the painters renderer to PNG (raster) and a
%   sibling SVG (vector).  When `path` ends in `.png`, an `.svg` is written
%   alongside it.  Both formats use the painters renderer for consistent
%   line-art rendering between raster and vector.
    set(fig, 'Renderer', 'painters', 'Color', 'w');
    drawnow;
    pause(0.2);
    drawnow;
    if endsWith(path, '.png')
        print(fig, path, '-dpng', sprintf('-r%d', dpi), '-painters');
        fprintf('  Saved: %s\n', path);
        svg_path = regexprep(path, '\.png$', '.svg');
        print(fig, svg_path, '-dsvg', '-painters');
        fprintf('  Saved: %s\n', svg_path);
    else
        exportgraphics(fig, path, 'Resolution', dpi, 'BackgroundColor', 'white');
        fprintf('  Saved: %s\n', path);
    end
    close(fig);
end
