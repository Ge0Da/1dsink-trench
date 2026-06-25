function s = computeArcLength(x, y)
%COMPUTEARCLENGTH Compute cumulative arc length along a path
%
%   s = computeArcLength(x, y)
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   VERSION 2.0 - Handles NaN values, validates metric coordinates
%
%   INPUTS:
%       x, y - Coordinate vectors (ordered, must be metric e.g. UTM meters)
%
%   OUTPUTS:
%       s    - Cumulative arc length (same size as x, y)
%              s(1) = 0, s(end) = total length
%              NaN propagates from invalid coordinates
%
%   NOTES:
%       - Assumes planar metric coordinates (e.g., UTM meters)
%       - If x,y are lon/lat degrees, results will be in degrees (meaningless)
%       - Single NaN in x or y makes s NaN from that point onward
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   NaN handling, coordinate validation

x = x(:);
y = y(:);

n = numel(x);
if n ~= numel(y)
    error('x and y must have same length');
end

s = zeros(n, 1);

if n < 2
    return;
end

% Check for lon/lat (warn if coordinates look like degrees)
if max(abs(x)) < 360 && max(abs(y)) < 90
    warning('computeArcLength: Coordinates appear to be lon/lat degrees. Results may be meaningless. Use projected coordinates (e.g., UTM meters).');
end

% Compute segment lengths
dx = diff(x);
dy = diff(y);
ds = sqrt(dx.^2 + dy.^2);

% Handle NaN: set segment to NaN if either endpoint is NaN
% This will propagate NaN through cumsum from that point
nan_mask = isnan(dx) | isnan(dy);
ds(nan_mask) = NaN;

% Cumulative sum
s = [0; cumsum(ds)];

% If there are NaNs, they will propagate. Warn user.
if any(nan_mask)
    nNaN = sum(nan_mask);
    warning('computeArcLength: %d NaN values in coordinates. Arc length is NaN from first NaN onward.', nNaN);
end

end
