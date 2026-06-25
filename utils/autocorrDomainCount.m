function [nMin, L_decorr, acf, lags_km] = autocorrDomainCount(s_km, signal, edgeMask)
%AUTOCORRDOMAINCOUNT  Data-driven minimum domain count from autocorrelation
%
%   [nMin, L_decorr, acf, lags_km] = autocorrDomainCount(s_km, signal, edgeMask)
%
%   Computes the decorrelation length of a segmentation signal (e.g. the
%   D_{2-20} curvature) from its autocorrelation function (ACF).  The
%   decorrelation length L_decorr is the lag at which the ACF first crosses
%   zero, representing the characteristic length scale of independent
%   morphological features.  The minimum plausible domain count is then:
%
%       nMin = floor( analysable_length / L_decorr )
%
%   This is fully data-driven — no arbitrary constants.
%
%   INPUTS
%       s_km      - Along-axis distance vector (km), uniformly sampled,
%                   monotonic increasing.  [N x 1] or [1 x N] double.
%       signal    - The segmentation signal (e.g. curvature in m/km^2),
%                   same length as s_km.  [N x 1] or [1 x N] double.
%       edgeMask  - Logical mask, true = valid interior point (outside
%                   edge-degradation zones).  Same length as s_km.
%
%   OUTPUTS
%       nMin      - Minimum plausible domain count (integer >= 1).
%       L_decorr  - Decorrelation length (km).
%       acf       - Normalised autocorrelation function for valid signal
%                   (positive lags only, acf(1) = 1).  [M x 1] double.
%       lags_km   - Lag values in km corresponding to acf.  [M x 1] double.
%
%   EDGE CASES
%       - ACF never crosses zero: uses the lag where ACF drops below 1/e
%         as L_decorr and issues a warning.
%       - L_decorr shorter than typical minDomainLength_km: returns nMin
%         based on L_decorr but warns that merging will reduce actual count.
%       - Signal all-NaN or all-zero in valid region: errors with message.
%
%   EXAMPLE
%       ld = load('output/1DSink_MasterFigure_results.mat');
%       [nMin, Ld, acf, lags] = autocorrDomainCount( ...
%           ld.bands.s_km, ld.bands.d2D_2_20_mpkm2, ld.bands.edgeMask);
%       fprintf('L_decorr = %.1f km, nMin = %d\n', Ld, nMin);
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0


%% ---- Input validation ------------------------------------------------
s_km     = s_km(:);
signal   = signal(:);
edgeMask = edgeMask(:);

n = numel(s_km);
assert(numel(signal) == n && numel(edgeMask) == n, ...
    'autocorrDomainCount:sizeMismatch', ...
    'All inputs must have the same length (got %d, %d, %d).', ...
    n, numel(signal), numel(edgeMask));

assert(islogical(edgeMask), ...
    'autocorrDomainCount:badMask', ...
    'edgeMask must be logical.');

%% ---- Extract valid (non-edge) signal ---------------------------------
sig_valid = signal(edgeMask);
s_valid   = s_km(edgeMask);

% Remove any internal NaN values (safety — should not occur in normal use)
nanMask   = isnan(sig_valid);
if any(nanMask)
    fprintf('  autocorrDomainCount: removing %d NaN values from valid region.\n', sum(nanMask));
    sig_valid = sig_valid(~nanMask);
    s_valid   = s_valid(~nanMask);
end

Nvalid = numel(sig_valid);

% Check for degenerate inputs
if Nvalid < 10
    error('autocorrDomainCount:insufficientData', ...
        'Fewer than 10 valid non-NaN points (got %d). Cannot compute ACF.', Nvalid);
end

if all(sig_valid == 0)
    error('autocorrDomainCount:zeroSignal', ...
        'Signal is identically zero in the valid region. Cannot compute ACF.');
end

if std(sig_valid) == 0
    error('autocorrDomainCount:constantSignal', ...
        'Signal is constant in the valid region. Cannot compute ACF.');
end

%% ---- Compute analysable length ----------------------------------------
analysable_length_km = s_valid(end) - s_valid(1);
ds_km = median(diff(s_valid));   % uniform spacing (should be ~0.1 km)

fprintf('  autocorrDomainCount: %d valid points, analysable length = %.1f km, spacing = %.3f km\n', ...
    Nvalid, analysable_length_km, ds_km);

%% ---- Compute normalised ACF via FFT (Wiener-Khinchin) -----------------
%  This avoids dependency on Signal Processing Toolbox (no xcorr needed).
%  Zero-mean the signal, zero-pad to avoid circular correlation artefacts,
%  then use: ACF = ifft(|fft(x)|^2) / var.

sig_zm = sig_valid - mean(sig_valid);

Nfft   = 2^nextpow2(2 * Nvalid - 1);   % zero-pad for linear correlation
S      = fft(sig_zm, Nfft);
acf_full = ifft(S .* conj(S), 'symmetric');

% Keep positive lags only (0 to Nvalid-1) and normalise so acf(1) = 1
acf     = acf_full(1:Nvalid) / acf_full(1);
lags_km = (0:Nvalid-1)' * ds_km;

%% ---- Find decorrelation length ----------------------------------------
%  L_decorr = first lag where ACF crosses zero.
%  If ACF never crosses zero, fall back to the 1/e threshold.

zeroIdx = find(acf(2:end) <= 0, 1, 'first');   % skip lag=0 (acf=1)

if ~isempty(zeroIdx)
    % Linear interpolation between the last positive and first non-positive
    idx_before = zeroIdx;       % acf(idx_before) > 0  (1-indexed, skipped lag 0)
    idx_after  = zeroIdx + 1;   % acf(idx_after)  <= 0
    % Interpolate to find exact zero-crossing lag
    a0 = acf(idx_before);
    a1 = acf(idx_after);
    frac = a0 / (a0 - a1);     % fraction of the interval to the zero crossing
    L_decorr = lags_km(idx_before) + frac * ds_km;
    fprintf('  autocorrDomainCount: ACF zero-crossing at lag = %.1f km\n', L_decorr);
else
    % ACF never crosses zero — use 1/e threshold
    eIdx = find(acf(2:end) <= exp(-1), 1, 'first');
    if ~isempty(eIdx)
        L_decorr = lags_km(eIdx + 1);
    else
        % ACF stays above 1/e for the entire profile — use half the profile
        L_decorr = analysable_length_km / 2;
    end
    warning('autocorrDomainCount:noZeroCrossing', ...
        ['ACF does not cross zero. Using 1/e fallback: L_decorr = %.1f km.\n' ...
         'This may indicate a strong monotonic trend in the signal.'], L_decorr);
end

%% ---- Compute minimum domain count ------------------------------------
if L_decorr <= 0
    error('autocorrDomainCount:badLdecorr', ...
        'Computed L_decorr = %.3f km (non-positive). Check input signal.', L_decorr);
end

nMin = floor(analysable_length_km / L_decorr);
nMin = max(nMin, 1);   % at least 1 domain

fprintf('  autocorrDomainCount: L_decorr = %.1f km, analysable = %.1f km → nMin = %d\n', ...
    L_decorr, analysable_length_km, nMin);

% Warn if L_decorr is very short (shorter than typical minDomainLength)
if L_decorr < 20
    warning('autocorrDomainCount:shortLdecorr', ...
        ['L_decorr = %.1f km is shorter than typical minDomainLength_km (20 km).\n' ...
         'Domain merging in segmentDomains will reduce the actual count below nMin = %d.'], ...
        L_decorr, nMin);
end

end
