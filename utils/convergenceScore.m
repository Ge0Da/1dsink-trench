function [composites, details] = convergenceScore(boundaryPositions_km, s_km, curvatureSignal, slopeSignal, roughnessSignal, windowHalfWidth_km)
%CONVERGENCESCORE  Continuous multi-proxy convergence score for domain boundaries
%
%   [composites, details] = convergenceScore(boundaryPositions_km, s_km, ...
%       curvatureSignal, slopeSignal, roughnessSignal, windowHalfWidth_km)
%
%   For each boundary position, computes contrast ratios (CR) for three
%   independent morphological proxies — curvature, slope, and roughness —
%   then combines them into a single continuous composite score via the
%   geometric mean.  No thresholding, no discrete scores.
%
%   The contrast ratio at a boundary measures how different the signal is
%   on either side: CR = max(mean_before, mean_after) / min(...).
%   CR >= 1 by definition; higher values indicate stronger transitions.
%   The geometric mean of the three CRs gives a single number that is
%   high only when ALL proxies show a transition at that location.
%
%   INPUTS
%       boundaryPositions_km  - [nB x 1] boundary positions along axis (km)
%       s_km                  - [N x 1] along-axis distance vector (km),
%                               uniformly sampled, monotonic increasing
%       curvatureSignal       - [N x 1] curvature profile (e.g.
%                               d2D_2_20_mpkm2).  Absolute values are
%                               taken internally.
%       slopeSignal           - [N x 1] slope magnitude profile (e.g.
%                               abs(dD_2_10_mpkm)).  Absolute values are
%                               taken internally if signed.
%       roughnessSignal       - [N x 1] roughness proxy profile (e.g.
%                               abs(D_1_10_m)).  Absolute values are taken
%                               internally if signed.
%       windowHalfWidth_km    - scalar, half-width of the before/after
%                               comparison window (km).  Default: 20.
%
%   OUTPUTS
%       composites - [nB x 1] geometric mean of the three CRs per boundary
%       details    - struct with diagnostic fields:
%           .curv_CR        [nB x 1] contrast ratios for curvature
%           .slope_CR       [nB x 1] contrast ratios for slope
%           .rough_CR       [nB x 1] contrast ratios for roughness
%           .composite      [nB x 1] same as composites (for convenience)
%           .pct_populated  [nB x 3] fraction of window populated
%                           (col 1 = curv, col 2 = slope, col 3 = rough)
%                           Each value = min(pct_before, pct_after).
%
%   ALGORITHM
%       For each boundary b_km:
%         1. Define before = [b_km - W, b_km], after = [b_km, b_km + W]
%         2. For each proxy, compute mean of |signal| in each window
%         3. CR = max(mean_before, mean_after) / min(mean_before, mean_after)
%         4. composite = (curv_CR * slope_CR * rough_CR)^(1/3)
%
%   EDGE CASES
%       - Window extends beyond profile: uses available data
%       - Proxy all-NaN in a window: CR = 1 (no contrast, conservative)
%       - Division by zero (mean = 0): CR = 1
%       - Zero boundaries: returns empty arrays
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0


%% ---- Input validation ------------------------------------------------
boundaryPositions_km = boundaryPositions_km(:);
s_km            = s_km(:);
curvatureSignal = curvatureSignal(:);
slopeSignal     = slopeSignal(:);
roughnessSignal = roughnessSignal(:);

N  = numel(s_km);
nB = numel(boundaryPositions_km);

assert(numel(curvatureSignal) == N, ...
    'convergenceScore:sizeMismatch', ...
    'curvatureSignal length (%d) must match s_km length (%d).', ...
    numel(curvatureSignal), N);
assert(numel(slopeSignal) == N, ...
    'convergenceScore:sizeMismatch', ...
    'slopeSignal length (%d) must match s_km length (%d).', ...
    numel(slopeSignal), N);
assert(numel(roughnessSignal) == N, ...
    'convergenceScore:sizeMismatch', ...
    'roughnessSignal length (%d) must match s_km length (%d).', ...
    numel(roughnessSignal), N);
assert(isscalar(windowHalfWidth_km) && windowHalfWidth_km > 0, ...
    'convergenceScore:badWindow', ...
    'windowHalfWidth_km must be a positive scalar.');

%% ---- Handle empty input ----------------------------------------------
if nB == 0
    composites = zeros(0, 1);
    details = struct('curv_CR', zeros(0,1), 'slope_CR', zeros(0,1), ...
        'rough_CR', zeros(0,1), 'composite', zeros(0,1), ...
        'pct_populated', zeros(0,3));
    return
end

%% ---- Take absolute values (signals may be signed) --------------------
curv_abs  = abs(curvatureSignal);
slope_abs = abs(slopeSignal);
rough_abs = abs(roughnessSignal);

%% ---- Compute contrast ratios for each boundary -----------------------
ds_km = median(diff(s_km));   % uniform spacing
expectedPoints = round(windowHalfWidth_km / ds_km);

curv_CR  = ones(nB, 1);
slope_CR = ones(nB, 1);
rough_CR = ones(nB, 1);
pct_pop  = zeros(nB, 3);   % min(pct_before, pct_after) per proxy

for i = 1:nB
    b = boundaryPositions_km(i);

    % Define window masks
    mask_before = (s_km >= b - windowHalfWidth_km) & (s_km <= b);
    mask_after  = (s_km >= b) & (s_km <= b + windowHalfWidth_km);

    % --- Curvature ---
    [curv_CR(i), pct_pop(i,1)] = computeCR(curv_abs, mask_before, mask_after, expectedPoints);

    % --- Slope ---
    [slope_CR(i), pct_pop(i,2)] = computeCR(slope_abs, mask_before, mask_after, expectedPoints);

    % --- Roughness ---
    [rough_CR(i), pct_pop(i,3)] = computeCR(rough_abs, mask_before, mask_after, expectedPoints);
end

%% ---- Composite: geometric mean of three CRs -------------------------
composites = (curv_CR .* slope_CR .* rough_CR) .^ (1/3);

%% ---- Pack details struct ---------------------------------------------
details.curv_CR       = curv_CR;
details.slope_CR      = slope_CR;
details.rough_CR      = rough_CR;
details.composite     = composites;
details.pct_populated = pct_pop;

%% ---- Print summary ---------------------------------------------------
fprintf('  convergenceScore: %d boundaries, composite CR range [%.3f, %.3f]\n', ...
    nB, min(composites), max(composites));

end


%% ========================================================================
%  LOCAL HELPER: compute contrast ratio for one proxy at one boundary
%  ========================================================================
function [CR, pct_min] = computeCR(proxyAbs, mask_before, mask_after, expectedPts)
%COMPUTECR  Contrast ratio between before and after windows
%
%   Returns CR >= 1 (ratio of larger mean to smaller mean).
%   Returns CR = 1 if data is insufficient or zero.
%   pct_min = min(pct_before, pct_after) for diagnostics.

    vals_before = proxyAbs(mask_before);
    vals_after  = proxyAbs(mask_after);

    % Remove NaN
    vals_before = vals_before(~isnan(vals_before));
    vals_after  = vals_after(~isnan(vals_after));

    n_b = numel(vals_before);
    n_a = numel(vals_after);

    % Population fractions
    if expectedPts > 0
        pct_b = n_b / expectedPts;
        pct_a = n_a / expectedPts;
    else
        pct_b = 0;
        pct_a = 0;
    end
    pct_min = min(pct_b, pct_a);

    % Compute means
    if isempty(vals_before) || isempty(vals_after)
        CR = 1;
        return
    end

    m_b = mean(vals_before);
    m_a = mean(vals_after);

    % Guard against division by zero
    if m_b == 0 && m_a == 0
        CR = 1;
    elseif min(m_b, m_a) == 0
        CR = Inf;
    else
        CR = max(m_b / m_a, m_a / m_b);
    end
end
