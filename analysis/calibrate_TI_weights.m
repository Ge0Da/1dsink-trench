% calibrate_TI_weights.m
% =========================================================================
%  TI weight calibration against particle-reactive heavy-metal geochemistry
% =========================================================================
%
%  Calibrates the composite Trapping Index (TI) weights against independent
%  sediment geochemistry. The calibration target is the mean normalized
%  heavy-metal content (Hg, As, Cd, Pb; n = 12 stations). Cadmium is retained
%  in the target although it correlates near zero with the TI, which shows
%  that particle-reactive morphological trapping is selective by metal
%  chemistry rather than granted for all heavy metals (Discussion).
%
%  Because three weights cannot be identified from twelve stations, the
%  grid-search optimum is reported only to quantify overfitting (its
%  leave-one-out cross-validation is far below its calibration-set fit). The
%  index is instead built from fixed, balanced a-priori weights (0.47, 0.40,
%  0.13) within the bootstrap envelope, validated against the geochemistry.
%
%  Steps:
%    1. Single-component diagnostics
%    2. Simplex grid search (weights sum to 1)
%    3. Leave-one-out cross-validation
%    4. Bootstrap confidence intervals
%    5. Save results (figures are built by paper_figures.m)
%
%  Output: output/TI_weight_calibration.mat
%
%  Author: David Fernandez-Blanco
% =========================================================================

if ~exist('tTotal', 'var'), clear; clc; end
projRoot = fileparts(fileparts(mfilename('fullpath')));
if isempty(projRoot), projRoot = fileparts(pwd); end
addpath(fullfile(projRoot, 'utils'));
cd(projRoot);

fprintf('\nTI weight calibration (four-element heavy-metal target)\n');
fprintf('%s\n', repmat('=', 1, 60));

%% 1. Load data ------------------------------------------------------------
ld = load('output/1DSink_MasterFigure_results.mat');
T  = readtable('data/kamchatka_samples.xlsx', 'VariableNamingRule', 'preserve');
colNames = T.Properties.VariableNames;

ps          = ld.projectedSamples;
sample_s_km = ps.s_km(:);
mr          = ld.masterResults;
metrics_s_km = mr.metrics.s_km(:);

% Component scores interpolated to sample positions (band-coarse manifold).
slope_at = interp1(metrics_s_km, mr.metrics.components.slope_norm, sample_s_km, 'linear', NaN);
rough_at = interp1(metrics_s_km, mr.metrics.components.rough_norm, sample_s_km, 'linear', NaN);
curv_at  = interp1(metrics_s_km, mr.metrics.components.curv_norm,  sample_s_km, 'linear', NaN);

% Operational TI at samples from the dense per-km curve.
TI_s_km = ld.bands.s_km(:);
TI_val  = mr.metrics.TI_full(:);
vTI     = ~isnan(TI_val) & ~isnan(TI_s_km);
TI_at   = interp1(TI_s_km(vTI), TI_val(vTI), sample_s_km, 'linear', NaN);

% Calibration target: mean normalized content of the metals in calib_config()
% (default Hg, As, Cd, Pb). Cadmium is retained although it correlates near
% zero with the TI (see Discussion and the Hg-As-Pb sensitivity test).
HM_target = targetMean(T, colNames);

validAll = ~isnan(slope_at) & ~isnan(rough_at) & ~isnan(curv_at) & ~isnan(HM_target);
S  = slope_at(validAll);
R  = rough_at(validAll);
C  = curv_at(validAll);
HM = HM_target(validAll);
TI_current = TI_at(validAll);
n = numel(HM);
fprintf('  %d stations with complete components + heavy-metal data\n\n', n);

%% 2. Single-component diagnostics ----------------------------------------
compNames = {'LowSlope', 'Roughness', 'Concavity'};
compVecs  = {S, R, C};
compRho   = nan(3, 1);
compP     = nan(3, 1);
for i = 1:3
    [compRho(i), compP(i)] = corr(compVecs{i}, HM, 'type', 'Spearman');
    fprintf('  %-10s rho = %+0.3f  p = %0.4f\n', compNames{i}, compRho(i), compP(i));
end
[rho_current, p_current] = corr(TI_current, HM, 'type', 'Spearman');
fprintf('  %-10s rho = %+0.3f  p = %0.4f   (operational TI, dense)\n\n', 'TI', rho_current, p_current);

%% 3. Simplex grid search --------------------------------------------------
step  = 0.01;
wvals = 0:step:1.0;
nw    = numel(wvals);
rhoGrid = NaN(nw, nw);
bestRho = -Inf; bestWs = [NaN NaN NaN];
for i = 1:nw
    ws = wvals(i);
    for j = 1:nw
        wr = wvals(j);
        wc = 1.0 - ws - wr;
        if wc < -1e-10, continue; end
        wc = max(0, wc);
        rho = corr(ws*S + wr*R + wc*C, HM, 'type', 'Spearman');
        rhoGrid(i, j) = rho;
        if rho > bestRho, bestRho = rho; bestWs = [ws, wr, wc]; end
    end
end
[~, bestP] = corr(bestWs(1)*S + bestWs(2)*R + bestWs(3)*C, HM, 'type', 'Spearman');

% Adopted operational weights (1.4 : 1.2 : 0.4 normalized), fixed a priori.
wCur = [1.4, 1.2, 0.4] / 3.0;
fprintf('  Grid-search optimum: w = [%.2f %.2f %.2f]  rho = %+0.3f (p = %0.4f)\n', bestWs, bestRho, bestP);
fprintf('  Operational weights: w = [%.3f %.3f %.3f]  rho = %+0.3f (dense)\n\n', wCur, rho_current);

%% 4. Leave-one-out cross-validation --------------------------------------
TI_loocv = nan(n, 1);
for k = 1:n
    msk = true(n, 1); msk(k) = false;
    Sk = S(msk); Rk = R(msk); Ck = C(msk); HMk = HM(msk);
    bestRho_k = -Inf; bestW_k = [1/3 1/3 1/3];
    for i = 1:nw
        ws = wvals(i);
        for j = 1:nw
            wr = wvals(j);
            wc = 1.0 - ws - wr;
            if wc < -1e-10, continue; end
            wc = max(0, wc);
            rho_k = corr(ws*Sk + wr*Rk + wc*Ck, HMk, 'type', 'Spearman');
            if rho_k > bestRho_k, bestRho_k = rho_k; bestW_k = [ws, wr, wc]; end
        end
    end
    TI_loocv(k) = bestW_k(1)*S(k) + bestW_k(2)*R(k) + bestW_k(3)*C(k);
end
[rho_loocv, p_loocv] = corr(TI_loocv, HM, 'type', 'Spearman');
fprintf('  LOOCV rho = %+0.3f   (optimum %+0.3f; %.0f%% shrinkage)\n\n', ...
    rho_loocv, bestRho, 100 * (bestRho - rho_loocv) / abs(bestRho));

%% 5. Bootstrap confidence intervals --------------------------------------
nBoot = 1000;
bootW   = nan(nBoot, 3);
bootRho = nan(nBoot, 1);
rng(42);
for b = 1:nBoot
    idx = randi(n, n, 1);
    Sb = S(idx); Rb = R(idx); Cb = C(idx); HMb = HM(idx);
    bestRho_b = -Inf; bestW_b = [1/3 1/3 1/3];
    for ws = 0:0.05:1.0
        for wr = 0:0.05:(1.0 - ws)
            wc = 1.0 - ws - wr;
            rho_b = corr(ws*Sb + wr*Rb + wc*Cb, HMb, 'type', 'Spearman');
            if rho_b > bestRho_b, bestRho_b = rho_b; bestW_b = [ws, wr, wc]; end
        end
    end
    bootW(b, :) = bestW_b;
    bootRho(b)  = bestRho_b;
end
ci_w   = prctile(bootW, [2.5, 97.5], 1);
med_w  = median(bootW, 1);
fprintf('  Bootstrap median w = [%.2f %.2f %.2f]\n', med_w);
fprintf('  Bootstrap 95%% CI w_slope [%.2f, %.2f]\n\n', ci_w(1,1), ci_w(2,1));

%% 6. Save results ---------------------------------------------------------
if ~exist('output', 'dir'), mkdir('output'); end
results = struct();
results.n_samples = n;
results.component_rho   = compRho;
results.component_p     = compP;
results.component_names = {compNames};
results.current_weights_normalised = wCur;
results.current_rho     = rho_current;
results.current_p       = p_current;
results.optimal_weights = bestWs;
results.optimal_rho     = bestRho;
results.optimal_p       = bestP;
results.rhoGrid         = rhoGrid;
results.wvals           = wvals;
results.loocv_rho       = rho_loocv;
results.loocv_p         = p_loocv;
results.loocv_TI        = TI_loocv;
results.shrinkage       = bestRho - rho_loocv;
results.bootstrap_weights        = bootW;
results.bootstrap_rho_all        = bootRho;
results.bootstrap_weights_median = med_w;
results.bootstrap_weights_ci95   = ci_w;
results.HM_norm_mean    = HM;
results.TI_current      = TI_current;
results.sample_s_km     = sample_s_km(validAll);
save('output/TI_weight_calibration.mat', 'results', '-v7.3');
fprintf('  Saved output/TI_weight_calibration.mat\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Calibration complete. Figures built by paper_figures.m.\n\n');

%% ------------------------------------------------------------------------
function v = getCol(Tbl, names, colname)
  if any(strcmp(names, colname))
    raw = Tbl.(colname);
    if iscell(raw)
      v = nan(numel(raw), 1);
      for jj = 1:numel(raw)
        if isnumeric(raw{jj}) && isscalar(raw{jj})
          v(jj) = double(raw{jj});
        elseif ischar(raw{jj}) || isstring(raw{jj})
          v(jj) = str2double(raw{jj});
        end
      end
    else
      v = double(raw);
    end
  else
    warning('calibrate_TI_weights: column "%s" not found.', colname);
    v = nan(height(Tbl), 1);
  end
end
