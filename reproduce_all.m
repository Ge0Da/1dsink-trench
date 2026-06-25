% reproduce_all.m
% Reproduces all results from a clean state.
% One command: reproduce_all
%
% Order:
%   1. Check requirements (MATLAB version, toolboxes, input data)
%   2. Run main pipeline (run_master_figure)
%   3. Run supplementary analyses (calibrate, geochem, sweep, stationarity, convergence)
%
% All output goes to output/. Previously generated files are overwritten.
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

clear; clc;
fprintf('\n');
fprintf('================================================================\n');
fprintf('  1DSink — Full Reproduction\n');
fprintf('================================================================\n\n');

%% ---- Check requirements -----------------------------------------------
fprintf('Checking requirements...\n');

% MATLAB version
v = ver('MATLAB');
vYear = str2double(v.Release(3:6));
if vYear < 2022
    error('reproduce_all: MATLAB R2022a or later required (found %s).', v.Release);
end
fprintf('  MATLAB %s — OK\n', v.Release);

% Required toolboxes
required = {'Mapping Toolbox', 'Signal Processing Toolbox', 'Statistics and Machine Learning Toolbox'};
installed = {ver().Name};
for i = 1:numel(required)
    if any(strcmp(installed, required{i}))
        fprintf('  %s — OK\n', required{i});
    else
        error('reproduce_all: missing toolbox: %s', required{i});
    end
end

% Input data
if ~exist('data/kamchatka_gmrt.tif', 'file')
    error(['reproduce_all: DEM not found (data/kamchatka_gmrt.tif).\n' ...
           'Download from https://www.gmrt.org/ (Kuril-Kamchatka region).']);
end
if ~exist('data/kamchatka_samples.xlsx', 'file')
    error('reproduce_all: sample table not found (data/kamchatka_samples.xlsx).');
end
fprintf('  Input data — OK\n');

% Ensure utils/ is on path
addpath('utils');
fprintf('  utils/ on path — OK\n');

fprintf('\nAll requirements met.\n\n');

%% ---- Step 1: Main pipeline --------------------------------------------
fprintf('================================================================\n');
fprintf('  Step 1/6: Main pipeline (run_master_figure)\n');
fprintf('================================================================\n\n');
tTotal = tic;
run_master_figure

%% ---- Step 2: TI weight calibration ------------------------------------
fprintf('\n================================================================\n');
fprintf('  Step 2/6: TI weight calibration\n');
fprintf('================================================================\n\n');
run('analysis/calibrate_TI_weights')

%% ---- Step 3: Geochemical cross-check ----------------------------------
fprintf('\n================================================================\n');
fprintf('  Step 3/6: Geochemical cross-check\n');
fprintf('================================================================\n\n');
run('analysis/geochem_crosscheck')

%% ---- Step 4: Parameter sweep -------------------------------------------
fprintf('\n================================================================\n');
fprintf('  Step 4/6: Parameter sweep (constrained)\n');
fprintf('================================================================\n\n');
run('analysis/parameter_sweep')

%% ---- Step 5: Stationarity test -----------------------------------------
fprintf('\n================================================================\n');
fprintf('  Step 5/6: Stationarity test\n');
fprintf('================================================================\n\n');
run('analysis/stationarity_test')

%% ---- Step 6: Convergence test ------------------------------------------
fprintf('\n================================================================\n');
fprintf('  Step 6/6: Convergence scoring verification\n');
fprintf('================================================================\n\n');
run('analysis/test_convergence')

%% ---- Done --------------------------------------------------------------
elapsed = toc(tTotal);
fprintf('\n================================================================\n');
fprintf('  REPRODUCTION COMPLETE — %.1f s (%.1f min)\n', elapsed, elapsed/60);
fprintf('  All outputs in: output/\n');
fprintf('================================================================\n');
