%% test_convergence.m
% =========================================================================
%  CONVEGENCE TEST — Stand alone
% =========================================================================
%
% Standalone test script for convergenceScore.m
%
% Loads the saved MAT file, extracts boundary positions and signal vectors,
% runs convergenceScore, and prints the continuous composite CR table.
% No discrete scores, no thresholds — just the raw contrast ratios and
% their geometric mean.
%
% Usage:
%   analysis/test_convergence
%   test_convergence
%
%
%  Author: David Fernández-Blanco
%  Date:   2026-02-26
% =========================================================================

% Only clear workspace when run standalone (not called from reproduce_all)
if ~exist('tTotal', 'var'), clear; clc; end

% Resolve project root (one level up from analysis/)
projRoot = fileparts(fileparts(mfilename('fullpath')));
if isempty(projRoot), projRoot = fileparts(pwd); end
addpath(fullfile(projRoot, 'utils'));
cd(projRoot);

% ─────────────────────────────────────────────────────────────────────────
%  1. Load data
% ─────────────────────────────────────────────────────────────────────────
% Try both filenames (project was renamed mid-development)
matFile = 'output/1DSink_MasterFigure_results.mat';
if ~exist(matFile, 'file')
    matFile = 'output/1DSink_MasterFigure_results.mat';
end
fprintf('Loading %s ...\n', matFile);
ld = load(matFile);

bands = ld.bands;
seg   = ld.precomputedSeg;
bpos  = seg.boundaries(:);

fprintf('Profile: %.0f km  (%d points)\n', bands.s_km(end) - bands.s_km(1), numel(bands.s_km));
fprintf('Boundaries: %d\n\n', numel(bpos));

% ─────────────────────────────────────────────────────────────────────────
%  2. Run convergenceScore
% ─────────────────────────────────────────────────────────────────────────
windowHW = 20;   % km
fprintf('Running convergenceScore (windowHalfWidth = %d km) ...\n\n', windowHW);

[composites, details] = convergenceScore( ...
    bpos, bands.s_km, ...
    bands.d2D_2_20_mpkm2, ...   % curvature
    bands.dD_2_10_mpkm, ...     % slope
    bands.D_1_10_m, ...         % roughness
    windowHW);

% ─────────────────────────────────────────────────────────────────────────
%  3. Print results table
% ─────────────────────────────────────────────────────────────────────────
fprintf('\n');
fprintf('  %-6s  %-10s  %-10s  %-10s  %-10s  %-12s  %-10s\n', ...
    'Bound', 'Pos_km', 'Curv_CR', 'Slope_CR', 'Rough_CR', 'Composite', 'PctPop_min');
fprintf('  %s\n', repmat('-', 1, 75));
for b = 1:numel(bpos)
    pctMin = min(details.pct_populated(b,:));
    fprintf('  B%-5d  %-10.1f  %-10.3f  %-10.3f  %-10.3f  %-12.3f  %-10.1f%%\n', ...
        b, bpos(b), details.curv_CR(b), details.slope_CR(b), ...
        details.rough_CR(b), details.composite(b), pctMin*100);
end
fprintf('  %s\n', repmat('-', 1, 75));

fprintf('\n  Summary:\n');
fprintf('    Composite CR: min=%.3f  median=%.3f  max=%.3f\n', ...
    min(composites), median(composites), max(composites));
fprintf('    Curvature CR: min=%.3f  median=%.3f  max=%.3f\n', ...
    min(details.curv_CR), median(details.curv_CR), max(details.curv_CR));
fprintf('    Slope CR:     min=%.3f  median=%.3f  max=%.3f\n', ...
    min(details.slope_CR), median(details.slope_CR), max(details.slope_CR));
fprintf('    Roughness CR: min=%.3f  median=%.3f  max=%.3f\n', ...
    min(details.rough_CR), median(details.rough_CR), max(details.rough_CR));

% ─────────────────────────────────────────────────────────────────────────
%  4. Compare with pre-computed convergence (if present in MAT)
% ─────────────────────────────────────────────────────────────────────────
if isfield(seg, 'convergence') && isfield(seg.convergence, 'composite')
    storedComp = seg.convergence.composite;
    maxDiff = max(abs(composites - storedComp));
    fprintf('\n  Cross-check vs stored convergence: max|diff| = %.6f', maxDiff);
    if maxDiff < 1e-10
        fprintf('  (MATCH)\n');
    else
        fprintf('  (MISMATCH — investigate)\n');
    end
end

fprintf('\ntest_convergence.m complete.\n');
