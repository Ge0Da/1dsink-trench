%% parameter_sweep.m
% =========================================================================
%  PARAMETER SWEEP: constrained segmentation sensitivity analysis
% =========================================================================
%  
% Constrained parameter sweep with data-driven domain count floor and
% physically-motivated ranking.
%
% Method:
%   - Data-driven domain-count floor from autocorrDomainCount, in place of a
%     fixed domain-number window.
%   - Persistence length capped at 20 km, consistent with the 20 km maximum
%     wavelength of the D_{2-20} segmentation band.
%   - Candidates ranked by number of robust boundaries (>= 0.5 stability),
%     with mean stability as the tiebreaker.
%   - Autocorrelation-derived constraint diagnostics are reported.
%
% Strategy: run the full grid WITHOUT stabilityTest (fast), filter to
% combos that produce >= nDomains_min domains, then re-run the top 10
% candidates WITH stabilityTest=true to obtain survival scores.
%
% Usage:
%   analysis/parameter_sweep
%   parameter_sweep
%
%
%  Author: David Fernández-Blanco
%  Date:   2026-03-09
% =========================================================================

% Only clear workspace when run standalone (not called from reproduce_all)
if ~exist('tTotal', 'var'), clear; clc; end
% Resolve project root (one level up from analysis/)
projRoot = fileparts(fileparts(mfilename('fullpath')));
if isempty(projRoot), projRoot = fileparts(pwd); end
addpath(fullfile(projRoot, 'utils'));
cd(projRoot);  % so output/ paths resolve correctly

% ─────────────────────────────────────────────────────────────────────────
%  1. Load pre-computed bands and config
% ─────────────────────────────────────────────────────────────────────────
fprintf('Loading output/1DSink_MasterFigure_results.mat ...\n');
ld = load('output/1DSink_MasterFigure_results.mat');

bands     = ld.bands;
cfg       = ld.cfg;
seg       = cfg.segmentation;        % shorthand

s_km      = bands.s_km;
curvature = bands.d2D_2_20_mpkm2;
edgeMask  = bands.edgeMask;

fprintf('Profile: %.0f km  (%d points)\n\n', s_km(end) - s_km(1), numel(s_km));

% ─────────────────────────────────────────────────────────────────────────
%  1b. Data-driven minimum domain count
% ─────────────────────────────────────────────────────────────────────────
fprintf('Computing data-driven domain count floor...\n');
[nDomains_min, L_decorr, acf_diag, lags_km_diag] = ...
    autocorrDomainCount(s_km, curvature, edgeMask);

fprintf('\n');
fprintf('  ┌─────────────────────────────────────────────────────────┐\n');
fprintf('  │  Autocorrelation decorrelation length: %7.1f km       │\n', L_decorr);
fprintf('  │  Data-driven minimum domain count:     %4d            │\n', nDomains_min);
fprintf('  └─────────────────────────────────────────────────────────┘\n');
fprintf('\n');

% ─────────────────────────────────────────────────────────────────────────
%  2. Define parameter grid
% ─────────────────────────────────────────────────────────────────────────
%  Persistence capped at 20 km.  The max wavelength in D_{2-20}
%  is 20 km, so persistence > 20 km is physically inconsistent with the
%  band — a zero-crossing that persists for longer than the longest
%  wavelength in the signal cannot be a band-internal feature.
sigmas_vals  = [1.5, 1.7, 2.0, 2.5];
persist_vals = [10,  20             ];
minDom_vals  = [20,  40,  60      ];
smooth_vals  = [10,  15,  20      ];
minStab_vals = [0,   0.3, 0.5     ];

nTotal = numel(sigmas_vals) * numel(persist_vals) * numel(minDom_vals) * ...
         numel(smooth_vals) * numel(minStab_vals);
fprintf('Parameter grid: %d combinations  (stabilityTest=false for speed)\n', nTotal);
fprintf('  persist_vals = [%s] (capped at 20 km, consistent with D_{2-20} band)\n\n', ...
    strjoin(arrayfun(@num2str, persist_vals, 'UniformOutput', false), ', '));

% ─────────────────────────────────────────────────────────────────────────
%  3. Allocate result arrays
% ─────────────────────────────────────────────────────────────────────────
res_sig  = zeros(nTotal, 1);
res_per  = zeros(nTotal, 1);
res_mnd  = zeros(nTotal, 1);
res_sw   = zeros(nTotal, 1);
res_ms   = zeros(nTotal, 1);
res_nD   = nan(nTotal, 1);
res_nB   = nan(nTotal, 1);
res_bpos = cell(nTotal, 1);   % variable-length; stored in cell

% ─────────────────────────────────────────────────────────────────────────
%  4. Main sweep loop
% ─────────────────────────────────────────────────────────────────────────
iC = 0;
tStart = tic;

for ds = sigmas_vals
  for pl = persist_vals
    for md = minDom_vals
      for sw = smooth_vals
        for ms = minStab_vals
          iC = iC + 1;
          res_sig(iC) = ds;
          res_per(iC) = pl;
          res_mnd(iC) = md;
          res_sw(iC)  = sw;
          res_ms(iC)  = ms;

          try
            [doms, bounds, ~] = segmentDomains(s_km, curvature, edgeMask, ...
              'deadbandSigmas',          ds, ...
              'persistenceLength_km',    pl, ...
              'minDomainLength_km',      md, ...
              'smoothWindow_km',         sw, ...
              'minStability',            ms, ...
              'deadbandMode',            seg.deadbandMode, ...
              'deadbandAmplitude',       seg.deadbandAmplitude_fixed, ...
              'deadbandFactor',          seg.deadbandFactor_percentile, ...
              'deadbandWindow_km',       seg.deadbandWindow_km, ...
              'confidenceWeights',       seg.confidenceWeights, ...
              'dwellScale_km',           seg.dwellScale_km, ...
              'localAmplitudeRadius_km', seg.localAmplitudeRadius_km, ...
              'maxBoundaries',           seg.maxBoundaries, ...
              'stabilityTest',           false, ...
              'verbose',                 false);

            res_nD(iC)   = numel(doms);
            res_nB(iC)   = numel(bounds);
            res_bpos{iC} = bounds(:)';
          catch ME
            warning('parameter_sweep: combo %d failed: %s', iC, ME.message);
          end

          if mod(iC, 100) == 0
            elapsed = toc(tStart);
            fprintf('  %d / %d  (%.1f s elapsed)\n', iC, nTotal, elapsed);
          end
        end
      end
    end
  end
end

nOK = sum(~isnan(res_nD));
fprintf('\nSweep complete: %d / %d combos succeeded  (%.1f s)\n\n', ...
  nOK, nTotal, toc(tStart));

% Quick summary of domain-count distribution
uD = unique(res_nD(~isnan(res_nD)));
fprintf('Domain-count distribution:\n');
for k = 1:numel(uD)
  fprintf('  %2d domains : %d combos\n', uD(k), sum(res_nD == uD(k)));
end
fprintf('\n');

% ─────────────────────────────────────────────────────────────────────────
%  5. Filter: data-driven domain count floor
% ─────────────────────────────────────────────────────────────────────────
%  Use the autocorrelation-derived
%  nDomains_min.  No upper bound — let the data decide.
% NOTE: The autocorrelation-derived nDomains_min measures signal frequency,
% not geological domain spacing. The persistence cap (<=20 km) already
% eliminates trivial plate-boundary-scale solutions (2-domain configs).
% Use a minimal floor of 3 to exclude degenerate cases.
nDomains_floor = 3;
ok = ~isnan(res_nD) & res_nD >= nDomains_floor;
fprintf('Combos with >= %d domains (persistence-constrained floor): %d\n', nDomains_floor, sum(ok));

idxOK = find(ok);
if isempty(idxOK)
  error('parameter_sweep: no usable candidates found. Inspect domain-count distribution above.');
end

% Sort candidates by number of domains (descending) for re-run selection
% — prefer more-resolved segmentations for the stability test
[~, sortOrder] = sort(res_nD(idxOK), 'descend');
idxOK = idxOK(sortOrder);

% ─────────────────────────────────────────────────────────────────────────
%  6. Re-run top 10 with stabilityTest=true
% ─────────────────────────────────────────────────────────────────────────
nRerun = min(10, numel(idxOK));
fprintf('\nRe-running top %d candidates with stabilityTest=true ...\n', nRerun);
fprintf('(Each sub-run: %d stability combinations)\n\n', seg.stabilityGrid^3);

rr_origIdx  = zeros(nRerun, 1);
rr_nD       = nan(nRerun, 1);
rr_nB       = nan(nRerun, 1);
rr_nRobust  = nan(nRerun, 1);
rr_meanStab = nan(nRerun, 1);
rr_bpos     = cell(nRerun, 1);
rr_doms     = cell(nRerun, 1);
rr_diag     = cell(nRerun, 1);

for k = 1:nRerun
  ii = idxOK(k);
  fprintf('  [%d/%d]  σ=%.1f  persist=%.0f km  minDomain=%.0f km  smooth=%.0f km  minStab=%.1f\n', ...
    k, nRerun, res_sig(ii), res_per(ii), res_mnd(ii), res_sw(ii), res_ms(ii));

  rr_origIdx(k) = ii;
  try
    [doms2, bounds2, diag2] = segmentDomains(s_km, curvature, edgeMask, ...
      'deadbandSigmas',          res_sig(ii), ...
      'persistenceLength_km',    res_per(ii), ...
      'minDomainLength_km',      res_mnd(ii), ...
      'smoothWindow_km',         res_sw(ii), ...
      'minStability',            res_ms(ii), ...
      'deadbandMode',            seg.deadbandMode, ...
      'deadbandAmplitude',       seg.deadbandAmplitude_fixed, ...
      'deadbandFactor',          seg.deadbandFactor_percentile, ...
      'deadbandWindow_km',       seg.deadbandWindow_km, ...
      'confidenceWeights',       seg.confidenceWeights, ...
      'dwellScale_km',           seg.dwellScale_km, ...
      'localAmplitudeRadius_km', seg.localAmplitudeRadius_km, ...
      'maxBoundaries',           seg.maxBoundaries, ...
      'stabilityTest',           true, ...
      'stabilityPerturbation',   seg.stabilityPerturbation, ...
      'stabilityGrid',           seg.stabilityGrid, ...
      'stabilitySnapRadius_km',  seg.stabilitySnapRadius_km, ...
      'verbose',                 false);

    stab             = diag2.boundaryStability;
    rr_nD(k)         = numel(doms2);
    rr_nB(k)         = numel(bounds2);
    rr_nRobust(k)    = sum(stab >= 0.5);
    rr_meanStab(k)   = mean(stab);
    rr_bpos{k}       = bounds2(:)';
    rr_doms{k}       = doms2;
    rr_diag{k}       = diag2;

    fprintf('    → nDomains=%d  nBoundaries=%d  nRobust(≥50%%)=%d  meanStab=%.3f\n', ...
      rr_nD(k), rr_nB(k), rr_nRobust(k), rr_meanStab(k));
  catch ME
    warning('parameter_sweep: re-run %d failed: %s', k, ME.message);
  end
end

% ─────────────────────────────────────────────────────────────────────────
%  7. Rank by nRobust (descending), tiebreaker: meanStab
% ─────────────────────────────────────────────────────────────────────────
%  Rank by number of robust boundaries (>= 0.5 stability),
%  not by mean stability.  This rewards more well-resolved boundaries
%  rather than fewer boundaries with trivially high mean stability.
%  Tiebreaker: higher mean stability among equal nRobust counts.
validRR  = ~isnan(rr_meanStab);
robMask  = validRR & (rr_nRobust >= 1);   % at least 1 robust boundary
fprintf('\nCandidates with >= %d domains AND >= 1 robust boundary: %d\n', ...
    nDomains_floor, sum(robMask));

if sum(robMask) == 0
  fprintf('  (relaxing — showing all successfully re-run combos)\n');
  robMask = validRR;
end

rankIdx = find(robMask);
% Sort by (-nRobust, -meanStab) so highest nRobust first, then highest meanStab
[~, sortRank] = sortrows([-rr_nRobust(rankIdx), -rr_meanStab(rankIdx)]);
rankIdx = rankIdx(sortRank);

% ─────────────────────────────────────────────────────────────────────────
%  8. Ranked table
% ─────────────────────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('─', 1, 90));
fprintf('RANKED PARAMETER SWEEP RESULTS  (sorted by nRobust desc, then meanStab desc)\n');
fprintf('%s\n', repmat('─', 1, 90));
fprintf('%-4s  %-6s  %-8s  %-7s  %-8s  %-8s  %-8s  %-8s  %-9s\n', ...
  'Rank', 'σ_MAD', 'Pers_km', 'mD_km', 'Smth_km', 'minStab', ...
  'nDomain', 'nRobust', 'meanStab');
fprintf('%s\n', repmat('─', 1, 90));

for k = 1:numel(rankIdx)
  r  = rankIdx(k);
  ii = rr_origIdx(r);
  fprintf('%-4d  %-6.1f  %-8.0f  %-7.0f  %-8.0f  %-8.1f  %-8d  %-8d  %-9.3f\n', ...
    k, res_sig(ii), res_per(ii), res_mnd(ii), res_sw(ii), res_ms(ii), ...
    rr_nD(r), rr_nRobust(r), rr_meanStab(r));
end
fprintf('%s\n\n', repmat('─', 1, 90));

% ─────────────────────────────────────────────────────────────────────────
%  9. Full domain listings for top 3 combos
% ─────────────────────────────────────────────────────────────────────────
topN = min(3, numel(rankIdx));
fprintf('=== FULL DOMAIN LISTINGS FOR TOP %d COMBOS ===\n\n', topN);

for k = 1:topN
  r  = rankIdx(k);
  ii = rr_origIdx(r);

  fprintf('┌── Rank %d ──────────────────────────────────────────────────────────\n', k);
  fprintf('│  σ_MAD=%.1f  persistenceLength=%.0f km  minDomainLength=%.0f km\n', ...
    res_sig(ii), res_per(ii), res_mnd(ii));
  fprintf('│  smoothWindow=%.0f km  minStability=%.1f\n', res_sw(ii), res_ms(ii));
  fprintf('│  nDomains=%d  nRobust=%d  meanStability=%.3f\n', ...
    rr_nD(r), rr_nRobust(r), rr_meanStab(r));
  fprintf('│\n');

  doms = rr_doms{r};
  dg   = rr_diag{r};

  if ~isempty(doms)
    stab = dg.boundaryStability;
    bpos = rr_bpos{r};

    fprintf('│  %-4s  %-10s  %-10s  %-9s  %-8s  %-10s  %-10s\n', ...
      'ID', 'Start_km', 'End_km', 'Len_km', 'Type', 'MeanCurv', 'Roughness');
    fprintf('│  %s\n', repmat('-', 1, 70));
    for d = 1:numel(doms)
      dm = doms(d);
      fprintf('│  %-4d  %-10.1f  %-10.1f  %-9.1f  %-8s  %-10.4f  %-10.4f\n', ...
        dm.id, dm.s_start, dm.s_end, dm.length_km, dm.type, ...
        dm.meanCurv, dm.meanRoughness);
    end

    fprintf('│\n│  Boundary stability:\n');
    for b = 1:numel(bpos)
      if b <= numel(stab)
        flag = '';
        if stab(b) >= 0.5
          flag = '  *** ROBUST';
        end
        fprintf('│    B%d: %6.1f km   survival=%.2f%s\n', b, bpos(b), stab(b), flag);
      end
    end

    % Check sink zone 650–900 km
    fprintf('│\n│  Sink-zone coverage (650–900 km expected major sink):\n');
    sinkFound = false;
    for d = 1:numel(doms)
      dm = doms(d);
      ov_s = max(dm.s_start, 650);
      ov_e = min(dm.s_end,   900);
      if ov_e > ov_s
        sinkFound = true;
        fprintf('│    Domain %d (%s, %.0f–%.0f km): %.0f km overlap\n', ...
          dm.id, dm.type, dm.s_start, dm.s_end, ov_e - ov_s);
      end
    end
    if ~sinkFound
      fprintf('│    (no domain overlaps 650–900 km — check profile range)\n');
    end
  else
    fprintf('│  (domain detail unavailable — re-run failed)\n');
  end
  fprintf('└────────────────────────────────────────────────────────────────────\n\n');
end

% ─────────────────────────────────────────────────────────────────────────
%  9b. Convergence scores for top-ranked configuration
% ─────────────────────────────────────────────────────────────────────────
%  Compute continuous multi-proxy convergence scores for the #1 config.
%  Uses curvature, slope, and roughness contrast ratios combined via
%  geometric mean — no thresholds, no discrete scores.
if ~isempty(rankIdx)
    topR  = rankIdx(1);
    topBounds = rr_bpos{topR};
    if ~isempty(topBounds)
        fprintf('=== CONVERGENCE SCORES FOR RANK-1 CONFIGURATION ===\n\n');
        convWindowHW = 20;   % km half-width
        [convComp, convDet] = convergenceScore( ...
            topBounds(:), s_km, ...
            bands.d2D_2_20_mpkm2, ...   % curvature
            bands.dD_2_10_mpkm, ...     % slope
            bands.D_1_10_m, ...         % roughness
            convWindowHW);

        fprintf('  %-8s  %-10s  %-10s  %-10s  %-12s  %-10s\n', ...
            'Bound', 'Pos_km', 'Curv_CR', 'Slope_CR', 'Rough_CR', 'Composite');
        fprintf('  %s\n', repmat('-', 1, 65));
        for b = 1:numel(topBounds)
            fprintf('  B%-6d  %-10.1f  %-10.3f  %-10.3f  %-10.3f  %-12.3f\n', ...
                b, topBounds(b), convDet.curv_CR(b), convDet.slope_CR(b), ...
                convDet.rough_CR(b), convDet.composite(b));
        end
        fprintf('\n  Composite CR: min=%.3f  median=%.3f  max=%.3f\n\n', ...
            min(convComp), median(convComp), max(convComp));
    end
end

% ─────────────────────────────────────────────────────────────────────────
% 10. Save results
% ─────────────────────────────────────────────────────────────────────────
outMat = 'output/parameter_sweep_results.mat';
save(outMat, ...
  'res_sig', 'res_per', 'res_mnd', 'res_sw', 'res_ms', ...
  'res_nD',  'res_nB',  'res_bpos', ...
  'rr_origIdx', 'rr_nD', 'rr_nB', 'rr_nRobust', 'rr_meanStab', ...
  'rr_bpos', 'rr_doms', 'rr_diag', 'rankIdx', ...
  'nDomains_min', 'L_decorr', 'acf_diag', 'lags_km_diag');
fprintf('Results saved to %s\n', outMat);
fprintf('\nparameter_sweep.m complete.\n');
