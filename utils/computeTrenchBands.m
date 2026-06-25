function out = computeTrenchBands(s_m, z_m, opts)
%COMPUTETRENCHBANDS Multi-band decomposition of along-trench depth profile
%
%   out = computeTrenchBands(s_m, z_m)
%   out = computeTrenchBands(s_m, z_m, 'Name', value, ...)
%
%   Decomposes along-trench depth into geologically meaningful wavelength bands:
%     - Plot A (S2S):      Multiple sub-10 km bands for roughness analysis
%     - Plot B (Tectonics): 10-50 km, 50-200 km, 200-1000 km, >1000 km
%
%   VERSION 3.0 - Enhanced bands based on signal propagation theory
%
%   S2S BAND DEFINITIONS:
%   ============================
%   <10 km:      Full high-pass residual (amplitude reference only)
%   1-10 km:     Original band-pass
%   2-10 km:     PRIMARY for slope - removes sub-km noise, stable derivatives
%   5-10 km:     PRIMARY for curvature - yields interpretable domains
%   2-20 km:     SECONDARY for curvature - domain boundary picking
%
%   RATIONALE (from S2S signal propagation theory):
%   ================================================
%   - Derivatives weight short wavelengths as k^n (n=1 for slope, n=2 for curvature)
%   - <10 km band includes sub-km content that dominates derivatives
%   - 2-10 km removes gridding/axis-jitter artifacts from slope
%   - 5-10 km yields curvature that behaves as "domains" not "static"
%   - 2-20 km provides stable zero-crossings for segmentation boundaries
%
%   TECTONIC BAND DEFINITIONS:
%   ==========================
%   10-50 km:    Structural segmentation, incoming-plate roughness imprints
%   50-200 km:   Margin-scale domains (accretion vs erosion)
%   200-1000 km: Plate/slab-scale forcing (flexure, slab geometry)
%   >1000 km:    Very long wavelength (regional tilt)
%
%   INPUTS:
%       s_m         - Arc-length along path (meters)
%       z_m         - Elevation at each point (meters, negative for bathymetry)
%
%   NAME-VALUE OPTIONS:
%       'resampleSpacing_m' - Uniform resampling interval (m), default: 100
%       'useDepth'          - If true, output depth (+down), default: true
%       'interpMethod'      - Interpolation: 'makima'|'pchip'|'linear', default: 'makima'
%       'edgeTrim_km'       - Trim edges to avoid filter artifacts (km), default: 20
%       'derivWindow_km'    - Light smoothing for derivatives (km), default: 2
%
%   OUTPUTS:
%       out.s_km            - Distance along axis (km), uniformly sampled
%       out.depth_raw_m     - Raw depth (uniformly resampled)
%
%       % S2S Band components (depth):
%       out.D_lt10_m        - High-pass <10 km (amplitude reference)
%       out.D_1_10_m        - Band-pass 1-10 km
%       out.D_2_10_m        - Band-pass 2-10 km (slope band)
%       out.D_5_10_m        - Band-pass 5-10 km (curvature band)
%       out.D_2_20_m        - Band-pass 2-20 km (segmentation band)
%
%       % Tectonic Band components:
%       out.D_10_50_m       - Band-pass 10-50 km (structural)
%       out.D_50_200_m      - Band-pass 50-200 km (margin domains)
%       out.D_200_1000_m    - Band-pass 200-1000 km (plate/slab scale)
%       out.D_gt1000_m      - Low-pass >1000 km (regional)
%
%       % S2S Derivatives:
%       out.dD_2_10_mpkm    - Slope from 2-10 km (PRIMARY)
%       out.d2D_5_10_mpkm2  - Curvature from 5-10 km (PRIMARY)
%       out.d2D_2_20_mpkm2  - Curvature from 2-20 km (segmentation)
%
%       % Additional S2S band:
%       out.dD_lt10_mpkm    - Slope from <10 km (full high-pass)
%       out.dD_1_10_mpkm    - Slope from 1-10 km
%       out.d2D_lt10_mpkm2  - Curvature from <10 km
%       out.d2D_1_10_mpkm2  - Curvature from 1-10 km
%
%       % Tectonic Derivatives:
%       out.dD_10_50_mpkm, out.d2D_10_50_mpkm2
%       out.dD_50_200_mpkm, out.d2D_50_200_mpkm2
%       out.dD_200_1000_mpkm, out.d2D_200_1000_mpkm2
%       out.dD_gt1000_mpkm, out.d2D_gt1000_mpkm2
%
%       out.lowpass         - Struct with all low-pass components
%       out.stats           - Statistics for each band
%       out.edgeMask        - Logical mask for valid (non-edge) points
%
%   REFERENCES:
%       - Romans et al. (2016) - S2S routing and accommodation
%       - Straub et al. (2020) - Signal propagation in sedimentary systems
%       - Wang & Bilek (2014) - Subducting relief and segmentation
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 3.0

%% Parse inputs
arguments
    s_m (:,1) double
    z_m (:,1) double
    opts.resampleSpacing_m (1,1) double {mustBePositive} = 100
    opts.useDepth (1,1) logical = true
    opts.interpMethod (1,:) char {mustBeMember(opts.interpMethod, {'makima', 'pchip', 'linear'})} = 'makima'
    opts.edgeTrim_km (1,1) double {mustBeNonnegative} = 20
    opts.derivWindow_km (1,1) double {mustBePositive} = 2
    opts.lambda_km struct = struct()   % Override wavelength cutoffs (see config_master_figure)
end

%% Wavelength cutoffs (km) - extended for v3.0
% Defaults are used unless overridden via opts.lambda_km from config
lambda = struct(...
    'L1', 1, ...      % 1 km
    'L2', 2, ...      % 2 km (slope lower bound)
    'L5', 5, ...      % 5 km (curvature lower bound)
    'L10', 10, ...    % 10 km (S2S upper bound)
    'L20', 20, ...    % 20 km (segmentation upper bound)
    'L50', 50, ...    % 50 km
    'L200', 200, ...  % 200 km
    'L1000', 1000);   % 1000 km

% Merge any user-supplied lambda values (allows partial override)
if ~isempty(fieldnames(opts.lambda_km))
    flds = fieldnames(opts.lambda_km);
    for iF = 1:numel(flds)
        lambda.(flds{iF}) = opts.lambda_km.(flds{iF});
    end
end

%% Step 1: Clean and enforce strictly increasing s
valid = isfinite(s_m) & isfinite(z_m);
s = s_m(valid);
z = z_m(valid);

[s, ia] = unique(s, 'sorted');  % Sort AND remove duplicates
z = z(ia);

n = numel(s);
if n < 100
    warning('computeTrenchBands: Need at least 100 valid points for multi-band analysis');
    out = createEmptyOutput();
    return;
end

%% Step 2: Convert to depth convention
if opts.useDepth
    D = -z;  % Depth positive downward
else
    D = z;
end

%% Step 3: Resample to uniform spacing
s_min = s(1);
s_max = s(end);
ds_m = opts.resampleSpacing_m;

s_uniform = (s_min:ds_m:s_max)';
n_uniform = numel(s_uniform);

if n_uniform < 100
    warning('computeTrenchBands: Uniform grid too short for multi-band analysis');
    out = createEmptyOutput();
    return;
end

% Interpolate depth to uniform grid (makima avoids pchip overshoot)
D_uniform = interp1(s, D, s_uniform, opts.interpMethod, 'extrap');

%% Step 4: Convert to km
s_km = s_uniform / 1000;
ds_km = ds_m / 1000;
length_km = s_km(end) - s_km(1);

%% Step 5: Compute low-pass filters at each wavelength cutoff
% Gaussian filter: σ = λ / (2π)
LP_1 = gaussLP(D_uniform, lambda.L1, ds_km);
LP_2 = gaussLP(D_uniform, lambda.L2, ds_km);      % NEW v3.0
LP_5 = gaussLP(D_uniform, lambda.L5, ds_km);      % NEW v3.0
LP_10 = gaussLP(D_uniform, lambda.L10, ds_km);
LP_20 = gaussLP(D_uniform, lambda.L20, ds_km);    % NEW v3.0
LP_50 = gaussLP(D_uniform, lambda.L50, ds_km);
LP_200 = gaussLP(D_uniform, lambda.L200, ds_km);
LP_1000 = gaussLP(D_uniform, lambda.L1000, ds_km);

%% Step 6: Form band-pass components by subtraction
% S2S bands (<10 km)
D_lt10 = D_uniform - LP_10;       % High-pass <10 km (amplitude reference)
D_1_10 = LP_1 - LP_10;            % Band-pass 1-10 km
D_2_10 = LP_2 - LP_10;            % Band-pass 2-10 km (slope) - NEW v3.0
D_5_10 = LP_5 - LP_10;            % Band-pass 5-10 km (curvature) - NEW v3.0
D_2_20 = LP_2 - LP_20;            % Band-pass 2-20 km (segmentation) - NEW v3.0

% Tectonic bands (>10 km)
D_10_50 = LP_10 - LP_50;          % Band-pass 10-50 km
D_50_200 = LP_50 - LP_200;        % Band-pass 50-200 km
D_200_1000 = LP_200 - LP_1000;    % Band-pass 200-1000 km
D_gt1000 = LP_1000;               % Low-pass >1000 km

%% Step 6b: Detrended depth for spectral slope analysis (NEW v3.1)
% For spectral slope (β) estimation, we need DETRENDED raw data, not band-passed.
% Computing β on band-passed data measures filter roll-off, not seafloor physics.
% Reference: Goff & Jordan (1988), JGR - PSD slope requires unfiltered data.
%
% Detrending options:
%   D_detrend_50  = raw - LP_50   (removes >50 km: flexure, tilt)
%   D_detrend_1000 = raw - LP_1000 (removes >1000 km: very broad)
%
% Default recommendation: LP_50 detrending preserves 2-50 km roughness
D_detrend_50 = D_uniform - LP_50;
D_detrend_1000 = D_uniform - LP_1000;

%% Step 7: Compute derivatives for each band
derivWin_samples = max(3, round(opts.derivWindow_km / ds_km));
if mod(derivWin_samples, 2) == 0
    derivWin_samples = derivWin_samples + 1;  % Ensure odd for centered filter
end

% S2S bands - various smoothing levels
[dD_lt10, d2D_lt10] = computeBandDerivatives(D_lt10, ds_km, derivWin_samples);
[dD_1_10, d2D_1_10] = computeBandDerivatives(D_1_10, ds_km, derivWin_samples);
[dD_2_10, d2D_2_10] = computeBandDerivatives(D_2_10, ds_km, derivWin_samples);  % NEW v3.0
[~, d2D_5_10] = computeBandDerivatives(D_5_10, ds_km, derivWin_samples);         % NEW v3.0

% κ_{2-20} for segmentation: MUST be smoother than κ_{5-10}
% Use 3x more smoothing for stable domain boundaries
segmentSmoothWin = derivWin_samples * 3;
if mod(segmentSmoothWin, 2) == 0
    segmentSmoothWin = segmentSmoothWin + 1;
end
[~, d2D_2_20] = computeBandDerivatives(D_2_20, ds_km, segmentSmoothWin);          % NEW v3.0

% Tectonic bands - standard gradient (already smooth from LP)
[dD_10_50, d2D_10_50] = computeBandDerivatives(D_10_50, ds_km, 3);
[dD_50_200, d2D_50_200] = computeBandDerivatives(D_50_200, ds_km, 3);
[dD_200_1000, d2D_200_1000] = computeBandDerivatives(D_200_1000, ds_km, 3);
[dD_gt1000, d2D_gt1000] = computeBandDerivatives(D_gt1000, ds_km, 3);

%% Step 8: Create edge mask (ignore edges for interpretation)
% FIX-12: Validate that edgeTrim_km covers the longest filter's influence.
% The Gaussian is truncated at 4σ; if edgeTrim_km < 4σ of the longest
% filter, edge artifacts will leak into the "valid" region.
maxSigma_km = lambda.L1000 / (2*pi);  % σ for the longest wavelength filter
maxHalfWidth_km = 4 * maxSigma_km;    % Gaussian truncation radius
if opts.edgeTrim_km < maxHalfWidth_km
    warning('computeTrenchBands:edgeTrimShort', ...
        ['edgeTrim_km (%.1f) is shorter than longest filter half-width (%.1f km, ' ...
         'from L1000=%.0f km). Edge artifacts may leak into valid region.'], ...
        opts.edgeTrim_km, maxHalfWidth_km, lambda.L1000);
end
edgeTrim_samples = round(opts.edgeTrim_km / ds_km);
edgeMask = true(n_uniform, 1);
edgeMask(1:edgeTrim_samples) = false;
edgeMask(end-edgeTrim_samples+1:end) = false;

%% Step 9: Compute statistics (using valid interior only)
stats = struct();
stats.lt10 = bandStats(D_lt10, dD_lt10, d2D_lt10, edgeMask);
stats.b1_10 = bandStats(D_1_10, dD_1_10, d2D_1_10, edgeMask);
stats.b2_10 = bandStats(D_2_10, dD_2_10, d2D_2_10, edgeMask);     % NEW v3.0
stats.b5_10 = bandStats(D_5_10, [], d2D_5_10, edgeMask);          % NEW v3.0
stats.b2_20 = bandStats(D_2_20, [], d2D_2_20, edgeMask);          % NEW v3.0
stats.b10_50 = bandStats(D_10_50, dD_10_50, d2D_10_50, edgeMask);
stats.b50_200 = bandStats(D_50_200, dD_50_200, d2D_50_200, edgeMask);
stats.b200_1000 = bandStats(D_200_1000, dD_200_1000, d2D_200_1000, edgeMask);
stats.gt1000 = bandStats(D_gt1000, dD_gt1000, d2D_gt1000, edgeMask);

%% Build output structure
out = struct();

% Coordinates
out.s_km = s_km;
out.s_m = s_uniform;
out.ds_km = ds_km;
out.length_km = length_km;
out.nPoints = n_uniform;

% Raw depth
out.depth_raw_m = D_uniform;

% S2S Band components (depth)
out.D_lt10_m = D_lt10;        % <10 km (amplitude reference)
out.D_1_10_m = D_1_10;        % 1-10 km
out.D_2_10_m = D_2_10;        % 2-10 km (slope band) - NEW v3.0
out.D_5_10_m = D_5_10;        % 5-10 km (curvature band) - NEW v3.0
out.D_2_20_m = D_2_20;        % 2-20 km (segmentation band) - NEW v3.0

% Detrended depth for spectral analysis (NEW v3.1)
out.D_detrend_50_m = D_detrend_50;    % raw - LP_50 (for β computation)
out.D_detrend_1000_m = D_detrend_1000; % raw - LP_1000 (broad detrending)

% Tectonic Band components (depth)
out.D_10_50_m = D_10_50;
out.D_50_200_m = D_50_200;
out.D_200_1000_m = D_200_1000;
out.D_gt1000_m = D_gt1000;

% S2S First derivatives (m/km)
out.dD_lt10_mpkm = dD_lt10;   % <10 km (reference)
out.dD_1_10_mpkm = dD_1_10;   % 1-10 km
out.dD_2_10_mpkm = dD_2_10;   % 2-10 km (PRIMARY) - NEW v3.0

% Tectonic First derivatives (m/km)
out.dD_10_50_mpkm = dD_10_50;
out.dD_50_200_mpkm = dD_50_200;
out.dD_200_1000_mpkm = dD_200_1000;
out.dD_gt1000_mpkm = dD_gt1000;

% S2S Second derivatives (m/km²)
out.d2D_lt10_mpkm2 = d2D_lt10;   % <10 km (reference)
out.d2D_1_10_mpkm2 = d2D_1_10;   % 1-10 km
out.d2D_5_10_mpkm2 = d2D_5_10;   % 5-10 km (PRIMARY) - NEW v3.0
out.d2D_2_20_mpkm2 = d2D_2_20;   % 2-20 km (segmentation) - NEW v3.0

% Tectonic Second derivatives (m/km²)
out.d2D_10_50_mpkm2 = d2D_10_50;
out.d2D_50_200_mpkm2 = d2D_50_200;
out.d2D_200_1000_mpkm2 = d2D_200_1000;
out.d2D_gt1000_mpkm2 = d2D_gt1000;

% Low-pass signals (for reference)
out.lowpass = struct('LP_1', LP_1, 'LP_2', LP_2, 'LP_5', LP_5, ...
    'LP_10', LP_10, 'LP_20', LP_20, 'LP_50', LP_50, ...
    'LP_200', LP_200, 'LP_1000', LP_1000);

% Edge mask and statistics
out.edgeMask = edgeMask;
out.stats = stats;

% Parameters used
out.params = struct(...
    'resampleSpacing_m', opts.resampleSpacing_m, ...
    'interpMethod', opts.interpMethod, ...
    'edgeTrim_km', opts.edgeTrim_km, ...
    'derivWindow_km', opts.derivWindow_km, ...
    'lambda', lambda);

end

%% ========================================================================
% HELPER FUNCTIONS
% =========================================================================

function y_lp = gaussLP(y, lambda_km, ds_km)
%GAUSSLP Gaussian low-pass filter with wavelength cutoff
    sigma_km = lambda_km / (2 * pi);
    sigma_samples = sigma_km / ds_km;

    if sigma_samples < 0.5
        y_lp = y;  % No filtering needed
        return;
    end

    % Create Gaussian kernel (truncate at 4*sigma)
    halfWidth = ceil(4 * sigma_samples);
    kernel = exp(-(-halfWidth:halfWidth).^2 / (2 * sigma_samples^2));
    kernel = kernel / sum(kernel);

    % Apply convolution
    y_lp = conv(y, kernel, 'same');

    % Fix edge effects with reflection padding
    nEdge = halfWidth;
    if nEdge > 0 && numel(y) > 2*nEdge
        % Reflect-pad and re-convolve edges
        y_padded = [flipud(y(2:nEdge+1)); y; flipud(y(end-nEdge:end-1))];
        y_conv = conv(y_padded, kernel, 'same');
        y_lp(1:nEdge) = y_conv(nEdge+1:2*nEdge);
        y_lp(end-nEdge+1:end) = y_conv(end-2*nEdge+1:end-nEdge);
    end
end

function [dD, d2D] = computeBandDerivatives(D_band, ds_km, smoothWin)
%COMPUTEBANDDERIVATIVES Compute first and second derivatives with optional smoothing

    % First derivative
    dD_raw = gradient(D_band, ds_km);  % m/km

    % Light smoothing for noisy bands
    if smoothWin > 1
        dD = smoothdata(dD_raw, 'movmean', smoothWin);
    else
        dD = dD_raw;
    end

    % Second derivative (from smoothed first derivative)
    d2D_raw = gradient(dD, ds_km);  % m/km²

    if smoothWin > 1
        d2D = smoothdata(d2D_raw, 'movmean', smoothWin);
    else
        d2D = d2D_raw;
    end
end

function s = bandStats(D, dD, d2D, mask)
%BANDSTATS Compute statistics for a band (interior only)
    D_int = D(mask);

    s = struct();
    s.D_range = [min(D_int), max(D_int)];
    s.D_std = std(D_int, 'omitnan');

    if ~isempty(dD)
        dD_int = dD(mask);
        s.dD_median = median(abs(dD_int), 'omitnan');
        s.dD_p99 = prctile(abs(dD_int), 99);
    else
        s.dD_median = NaN;
        s.dD_p99 = NaN;
    end

    if ~isempty(d2D)
        d2D_int = d2D(mask);
        s.d2D_median = median(abs(d2D_int), 'omitnan');
        s.d2D_p99 = prctile(abs(d2D_int), 99);
    else
        s.d2D_median = NaN;
        s.d2D_p99 = NaN;
    end
end

function out = createEmptyOutput()
%CREATEEMPTYOUTPUT Return empty output structure
    out = struct();
    out.s_km = [];
    out.s_m = [];
    out.ds_km = NaN;
    out.length_km = NaN;
    out.nPoints = 0;
    out.depth_raw_m = [];

    % S2S bands
    out.D_lt10_m = [];
    out.D_1_10_m = [];
    out.D_2_10_m = [];
    out.D_5_10_m = [];
    out.D_2_20_m = [];

    % Detrended for spectral analysis (NEW v3.1)
    out.D_detrend_50_m = [];
    out.D_detrend_1000_m = [];

    % Tectonic bands
    out.D_10_50_m = [];
    out.D_50_200_m = [];
    out.D_200_1000_m = [];
    out.D_gt1000_m = [];

    % S2S derivatives
    out.dD_lt10_mpkm = [];
    out.dD_1_10_mpkm = [];
    out.dD_2_10_mpkm = [];
    out.d2D_lt10_mpkm2 = [];
    out.d2D_1_10_mpkm2 = [];
    out.d2D_5_10_mpkm2 = [];
    out.d2D_2_20_mpkm2 = [];

    % Tectonic derivatives
    out.dD_10_50_mpkm = [];
    out.dD_50_200_mpkm = [];
    out.dD_200_1000_mpkm = [];
    out.dD_gt1000_mpkm = [];
    out.d2D_10_50_mpkm2 = [];
    out.d2D_50_200_mpkm2 = [];
    out.d2D_200_1000_mpkm2 = [];
    out.d2D_gt1000_mpkm2 = [];

    out.lowpass = struct();
    out.edgeMask = [];
    out.stats = struct();
    out.params = struct();
end
