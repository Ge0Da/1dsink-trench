function metrics = computeSinkMetrics(bands, opts)
%COMPUTESINKMETRICS Compute sink strength metrics and Trapping Index
%
%   metrics = computeSinkMetrics(bands)
%   metrics = computeSinkMetrics(bands, 'windowSize_km', 50, ...)
%
%   STATIONARITY GUARANTEE
%   ==============================
%   The Trapping Index (TI) at every window depends ONLY on the data
%   inside that window (and its immediate neighbours for slope/curvature).
%   TI normalization uses ABSOLUTE PHYSICAL THRESHOLDS by default, not
%   global percentiles.  Trimming the profile does not change any TI value.
%
%   Percentile normalization is available via normMode='percentile'
%   but is documented as NON-STATIONARY: the same window can get different
%   TI scores depending on what other windows are in the profile.
%
%   SPECTRAL SLOPE (beta) METHODOLOGY
%   ==========================================
%   beta is the negative slope of PSD in log-log space: PSD(f) ~ f^(-beta)
%   Computed on DETRENDED raw depth (raw - LP_50), not band-passed data.
%
%   TRAPPING INDEX
%   ==============
%   TI = (w_slope * LowSlope + w_rough * Roughness + w_curv * Concave) / sum(w)
%
%   Default normalization (absolute thresholds):
%     - slope:  0 m/km -> score 1 (flat = trapping),  slopeHi m/km -> score 0
%     - rough:  0 m    -> score 0 (smooth = no trap),  roughHi m   -> score 1
%     - curv:   curvLo m/km^2 -> score 0 (convex),     curvHi m/km^2 -> score 1 (concave)
%
%   These thresholds have physical meaning: a slope of 15 m/km is steep
%   for a trench floor regardless of what other slopes exist in the profile.
%
%   INPUTS
%       bands   - Output from computeTrenchBands()
%
%   OPTIONS
%       windowSize_km    (50)    Moving window size (km)
%       windowOverlap    (0.8)   Overlap fraction
%       w_slope          (1.0)   TI weight for low-slope
%       w_rough          (1.2)   TI weight for roughness
%       w_curv           (0.8)   TI weight for curvature
%       normMode         ('absolute')  'absolute' (stationary) or 'percentile'
%       slopeLo_mpkm     (0.5)   Slope below this -> score 1.0
%       slopeHi_mpkm     (15.0)  Slope above this -> score 0.0
%       roughLo_m        (2.0)   RMS below this -> score 0.0
%       roughHi_m        (80.0)  RMS above this -> score 1.0
%       curvLo_mpkm2     (-0.01) Concavity below this -> score 0.0  (convex)
%       curvHi_mpkm2     (0.05)  Concavity above this -> score 1.0  (concave)
%       pctileLow        (5)     [percentile mode only]
%       pctileHigh       (95)    [percentile mode only]
%       betaDetrend      ('lp50')
%       betaFitBand      ([0.1, 0.5]) cyc/km
%
%   Part of 1DSink Along-Trench Analysis Pipeline

arguments
    bands struct
    opts.windowSize_km  (1,1) double {mustBePositive}  = 50
    opts.windowOverlap  (1,1) double {mustBeInRange(opts.windowOverlap,0,0.99)} = 0.8
    opts.w_slope        (1,1) double {mustBePositive}  = 1.0
    opts.w_rough        (1,1) double {mustBePositive}  = 1.2
    opts.w_curv         (1,1) double {mustBePositive}  = 0.8
    opts.normMode       (1,:) char {mustBeMember(opts.normMode,{'absolute','percentile'})} = 'absolute'
    opts.slopeLo_mpkm   (1,1) double {mustBeNonnegative} = 0.5
    opts.slopeHi_mpkm   (1,1) double {mustBePositive}    = 15.0
    opts.roughLo_m      (1,1) double {mustBeNonnegative}  = 2.0
    opts.roughHi_m      (1,1) double {mustBePositive}     = 80.0
    opts.curvLo_mpkm2   (1,1) double = -0.01
    opts.curvHi_mpkm2   (1,1) double = 0.05
    opts.pctileLow      (1,1) double {mustBeInRange(opts.pctileLow,0,50)}   = 5
    opts.pctileHigh     (1,1) double {mustBeInRange(opts.pctileHigh,50,100)} = 95
    opts.betaDetrend    (1,:) char {mustBeMember(opts.betaDetrend,{'lp50','lp1000'})} = 'lp50'
    opts.betaFitBand    (1,2) double {mustBePositive} = [0.1, 0.5]
end

%% ---- check fields -------------------------------------------------------
hasV3 = isfield(bands,'D_2_10_m') && isfield(bands,'dD_2_10_mpkm') ...
      && isfield(bands,'d2D_5_10_mpkm2');
if ~hasV3
    error('computeSinkMetrics:missingFields', ...
          'bands missing v3.0 fields.  Run updated computeTrenchBands.');
end

%% ---- extract data -------------------------------------------------------
s_km     = bands.s_km;
ds_km    = bands.ds_km;
edgeMask = bands.edgeMask;

edgeValidStart = find(edgeMask,1,'first');
edgeValidEnd   = find(edgeMask,1,'last');
if isempty(edgeValidStart), edgeValidStart = 1; end
if isempty(edgeValidEnd),   edgeValidEnd = numel(edgeMask); end

D_2_10     = bands.D_2_10_m;
slope_2_10 = bands.dD_2_10_mpkm;
curv_5_10  = bands.d2D_5_10_mpkm2;

%% ---- detrended signal for beta ------------------------------------------
switch opts.betaDetrend
    case 'lp50'
        if isfield(bands,'D_detrend_50_m')
            D_spectral = bands.D_detrend_50_m;
        elseif isfield(bands,'lowpass') && isfield(bands.lowpass,'LP_50')
            D_spectral = bands.depth_raw_m - bands.lowpass.LP_50;
        else
            error('computeSinkMetrics:noDetrend','Missing D_detrend_50_m and LP_50');
        end
    case 'lp1000'
        if isfield(bands,'D_detrend_1000_m')
            D_spectral = bands.D_detrend_1000_m;
        elseif isfield(bands,'lowpass') && isfield(bands.lowpass,'LP_1000')
            D_spectral = bands.depth_raw_m - bands.lowpass.LP_1000;
        else
            error('computeSinkMetrics:noDetrend','Missing D_detrend_1000_m and LP_1000');
        end
end

%% ---- window parameters --------------------------------------------------
windowSize_samples = round(opts.windowSize_km / ds_km);
stepSize_samples   = max(1, round(windowSize_samples * (1 - opts.windowOverlap)));
n    = length(s_km);
nWin = max(1, floor((n - windowSize_samples) / stepSize_samples) + 1);

win_centers = zeros(nWin,1);
rms_vals    = zeros(nWin,1);
mad_vals    = zeros(nWin,1);
beta_vals   = zeros(nWin,1);
slope_mag   = zeros(nWin,1);
curv_mean   = zeros(nWin,1);

%% ---- Welch parameters for beta ------------------------------------------
fmin_beta = opts.betaFitBand(1);
fmax_beta = opts.betaFitBand(2);
targetMaxWL = 1 / fmin_beta;
welchSegLen = min(round(1.5 * targetMaxWL / ds_km), windowSize_samples);
welchOvlp   = round(welchSegLen / 2);
welchNfft   = max(256, 2^nextpow2(welchSegLen));

%% ---- compute metrics in windows -----------------------------------------
for i = 1:nWin
    sI = (i-1)*stepSize_samples + 1;
    eI = min(sI + windowSize_samples - 1, n);

    if sI < edgeValidStart || eI > edgeValidEnd
        win_centers(i) = mean(s_km(sI:eI));
        rms_vals(i)  = NaN; mad_vals(i)  = NaN;
        beta_vals(i) = NaN; slope_mag(i) = NaN; curv_mean(i) = NaN;
        continue
    end

    win_centers(i) = mean(s_km(sI:eI));

    D_win    = D_2_10(sI:eI);
    slp_win  = slope_2_10(sI:eI);
    crv_win  = curv_5_10(sI:eI);
    Dsp_win  = D_spectral(sI:eI);

    Dv   = D_win(~isnan(D_win));
    slpv = slp_win(~isnan(slp_win));
    crvv = crv_win(~isnan(crv_win));
    Dspv = Dsp_win(~isnan(Dsp_win));

    if length(Dv) > 10
        rms_vals(i)  = sqrt(mean(Dv.^2));
        mad_vals(i)  = median(abs(Dv - median(Dv)));
        slope_mag(i) = median(abs(slpv));
        curv_mean(i) = mean(crvv);

        if length(Dspv) > 20
            try
                Dspd = Dspv - mean(Dspv);
                [psd,freq] = pwelch(Dspd, welchSegLen, welchOvlp, welchNfft, 1/ds_km);
                vIdx = freq >= fmin_beta & freq <= fmax_beta & psd > 0;
                if sum(vIdx) >= 4
                    p = polyfit(log10(freq(vIdx)), log10(psd(vIdx)), 1);
                    beta_vals(i) = -p(1);
                else
                    beta_vals(i) = NaN;
                end
            catch
                beta_vals(i) = NaN;
            end
        else
            beta_vals(i) = NaN;
        end
    else
        rms_vals(i) = NaN; mad_vals(i) = NaN;
        beta_vals(i) = NaN; slope_mag(i) = NaN; curv_mean(i) = NaN;
    end
end

%% ========================================================================
%  TRAPPING INDEX — SPATIALLY STATIONARY NORMALIZATION
%  ====================================================================
%
%  The core stationarity requirement: the TI score at any window must
%  depend ONLY on the data inside that window, not on the rest of the
%  profile.  This is achieved by normalizing against ABSOLUTE physical
%  thresholds with fixed, interpretable scales.
%
%  Why absolute?  A trench-floor slope of 15 m/km is steep regardless of
%  whether the adjacent 500 km of trench contains steeper or flatter
%  segments.  Global percentiles break this: adding a very steep canyon
%  to the profile would lower the TI at all previously-flat windows,
%  which is physically nonsensical.

switch opts.normMode
    case 'absolute'
        % ---- Slope: LOW slope = HIGH score (inverted) ----
        %  slopeLo -> score 1.0 (flat, promotes trapping)
        %  slopeHi -> score 0.0 (steep, promotes bypass)
        sR = opts.slopeHi_mpkm - opts.slopeLo_mpkm;
        if sR > 0
            slope_norm = 1 - (slope_mag - opts.slopeLo_mpkm) / sR;
        else
            slope_norm = 0.5 * ones(size(slope_mag));
        end

        % ---- Roughness: HIGH roughness = HIGH score ----
        %  roughLo -> score 0.0 (smooth)
        %  roughHi -> score 1.0 (rough, more mini-basins)
        rR = opts.roughHi_m - opts.roughLo_m;
        if rR > 0
            rough_norm = (rms_vals - opts.roughLo_m) / rR;
        else
            rough_norm = 0.5 * ones(size(rms_vals));
        end

        % ---- Concavity: positive concavity = HIGH score ----
        %  (sign convention: curv_concave positive when concave-up)
        curv_concave = -curv_mean;
        cR = opts.curvHi_mpkm2 - opts.curvLo_mpkm2;
        if cR > 0
            curv_norm = (curv_concave - opts.curvLo_mpkm2) / cR;
        else
            curv_norm = 0.5 * ones(size(curv_concave));
        end

    case 'percentile'
        % Percentile normalisation: non-stationary alternative to 'absolute'.
        % WARNING: TI values depend on which windows are in the profile.
        % Trimming the profile will change TI at untrimmed windows.
        warning('computeSinkMetrics:nonStationary', ...
            ['normMode=''percentile'' is NON-STATIONARY: TI values depend ' ...
             'on the full profile extent. Use ''absolute'' for trim-invariant results.']);

        pLo = opts.pctileLow;
        pHi = opts.pctileHigh;

        sv = slope_mag(~isnan(slope_mag));
        if ~isempty(sv)
            p5s = prctile(sv, pLo); p95s = prctile(sv, pHi);
            if p95s > p5s, slope_norm = 1 - (slope_mag - p5s)/(p95s - p5s);
            else,          slope_norm = 0.5 * ones(size(slope_mag)); end
        else
            slope_norm = NaN(size(slope_mag));
        end

        rv = rms_vals(~isnan(rms_vals));
        if ~isempty(rv)
            p5r = prctile(rv, pLo); p95r = prctile(rv, pHi);
            if p95r > p5r, rough_norm = (rms_vals - p5r)/(p95r - p5r);
            else,          rough_norm = 0.5 * ones(size(rms_vals)); end
        else
            rough_norm = NaN(size(rms_vals));
        end

        curv_concave = -curv_mean;
        cv = curv_concave(~isnan(curv_concave));
        if ~isempty(cv)
            p5c = prctile(cv, pLo); p95c = prctile(cv, pHi);
            if p95c > p5c, curv_norm = (curv_concave - p5c)/(p95c - p5c);
            else,          curv_norm = 0.5 * ones(size(curv_concave)); end
        else
            curv_norm = NaN(size(curv_concave));
        end
end

% Clip all components to [0, 1]
slope_norm = max(0, min(1, slope_norm));
rough_norm = max(0, min(1, rough_norm));
curv_norm  = max(0, min(1, curv_norm));

% Weighted combination
w_sum = opts.w_slope + opts.w_rough + opts.w_curv;
TI = (opts.w_slope * slope_norm + opts.w_rough * rough_norm + opts.w_curv * curv_norm) / w_sum;

%% ---- Interpolate TI to full resolution ----------------------------------
TI_full = interp1(win_centers, TI, s_km, 'pchip', NaN);

%% ========================================================================
%  HEATMAP NORMALIZATION  (display-only, does NOT affect analysis)
%  ====================================================================
%  NOTE: normalizeRow uses global min/max of each metric vector.
%  This is acceptable because:
%    1. It affects ONLY the color scale in heatmap visualizations
%    2. Raw metric values (metrics.rms, .mad, .beta) are unaffected
%    3. TI computation uses absolute thresholds (above), not heatmap colors
%  If you need trim-invariant heatmap colors, pass fixed color ranges
%  to the plotting function instead of relying on these normalized values.

rms_norm  = normalizeRowDisplay(rms_vals);
mad_norm  = normalizeRowDisplay(mad_vals);
beta_norm = normalizeRowDisplay(beta_vals);

%% ---- build output -------------------------------------------------------
metrics = struct();

metrics.s_km      = win_centers;
metrics.s_km_full = s_km;
metrics.nWindows  = nWin;

metrics.rms       = rms_vals;
metrics.mad       = mad_vals;
metrics.beta      = beta_vals;
metrics.slope_mag = slope_mag;
metrics.curv_mean = curv_mean;

metrics.rms_norm  = rms_norm;
metrics.mad_norm  = mad_norm;
metrics.beta_norm = beta_norm;

metrics.TI        = TI;
metrics.TI_full   = TI_full;

metrics.components = struct('slope_norm', slope_norm, ...
    'rough_norm', rough_norm, 'curv_norm', curv_norm);

metrics.params = struct( ...
    'windowSize_km',   opts.windowSize_km, ...
    'windowOverlap',   opts.windowOverlap, ...
    'w_slope',         opts.w_slope, ...
    'w_rough',         opts.w_rough, ...
    'w_curv',          opts.w_curv, ...
    'normMode',        opts.normMode, ...
    'slopeLo_mpkm',    opts.slopeLo_mpkm, ...
    'slopeHi_mpkm',    opts.slopeHi_mpkm, ...
    'roughLo_m',       opts.roughLo_m, ...
    'roughHi_m',       opts.roughHi_m, ...
    'curvLo_mpkm2',    opts.curvLo_mpkm2, ...
    'curvHi_mpkm2',    opts.curvHi_mpkm2, ...
    'betaDetrend',     opts.betaDetrend, ...
    'betaFitBand_cyc_km', opts.betaFitBand, ...
    'betaFitBand_wavelength_km', [1/opts.betaFitBand(2), 1/opts.betaFitBand(1)], ...
    'betaMethod',      'Welch PSD on detrended raw depth (Goff & Jordan 1988)');

end

%% ========================================================================
%  HELPER
%% ========================================================================

function norm = normalizeRowDisplay(vals)
%NORMALIZEROWDISPLAY  Min-max normalization for DISPLAY ONLY.
%   Uses global min/max of the input vector.  This is NOT spatially
%   stationary: trimming the profile changes the color scale.
%   Use absolute color ranges in the plotting function if you need
%   trim-invariant heatmap colors.
    vals = vals(:);
    v = vals(~isnan(vals));
    if ~isempty(v) && max(v) > min(v)
        norm = (vals - min(v)) / (max(v) - min(v));
        norm = max(0, min(1, norm));
    else
        norm = 0.5 * ones(size(vals));
    end
end
