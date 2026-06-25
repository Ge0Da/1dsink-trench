function [xo, yo] = orderPointsByConnectivity(x, y, cellsize)
%ORDERPOINTSBYCONNECTIVITY Order unordered points into continuous chain
%
%   [xo, yo] = orderPointsByConnectivity(x, y, cellsize)
%
%   Part of S2S ModuleA - Source-to-Sink Geomorphometric Analysis
%
%   Uses PCA to find principal axis, then sorts along it.
%   Handles curved/diagonal trenches that simple X/Y sorting fails on.
%   Preserves ALL original points (no interpolation).
%
%   WHY PCA INSTEAD OF SIMPLE X/Y SORTING:
%   ======================================
%   Simple X/Y sorting assumes the feature aligns with coordinate axes.
%   A trench running NE-SW will have overlapping X values, causing sort()
%   to create massive jumps between consecutive points. 
%   
%   Example failure: 549,745 points with "length" of 141,000,000 km 
%   (should be ~2,200 km) because sorting by X created huge Y jumps.
%
%   PCA finds the actual principal direction of the point cloud and
%   sorts along that axis, regardless of orientation.
%
%   SCIENTIFIC NOTE:
%   This function preserves ALL input points on the divide - no artificial
%   interpolation or subsampling that would create points not on the actual
%   drainage divide.
%
%   INPUTS:
%       x, y     - Coordinate vectors (unordered)
%       cellsize - Grid cell size (for duplicate removal threshold)
%
%   OUTPUTS:
%       xo, yo   - Ordered coordinate vectors (column vectors)
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   PCA-based sorting

x = x(:);
y = y(:);
n = numel(x);

if n < 2
    xo = x;
    yo = y;
    return;
end

% =========================================================================
% PCA-BASED SORTING
% =========================================================================
% Strategy: Find the principal axis of the point cloud via PCA,
% then sort by projection onto that axis. This handles curved features
% better than simple X/Y sorting.

% Center the data
xc = x - mean(x);
yc = y - mean(y);

% Compute covariance matrix
C = [xc, yc]' * [xc, yc] / n;

% Find eigenvectors/eigenvalues
[V, D] = eig(C);

% Principal axis is eigenvector with largest eigenvalue
[~, maxIdx] = max(diag(D));
principal_axis = V(:, maxIdx);

% Project all points onto principal axis
projection = xc * principal_axis(1) + yc * principal_axis(2);

% Sort by projection value
[~, sortIdx] = sort(projection);

xo = x(sortIdx);
yo = y(sortIdx);

% =========================================================================
% DUPLICATE REMOVAL
% =========================================================================
% Remove near-duplicates (same cell) while preserving order
% This is physically meaningful - duplicate coordinates add no information

if n > 1
    dx = diff(xo);
    dy = diff(yo);
    ds = hypot(dx, dy);
    
    % Keep points separated by at least a small fraction of cellsize
    minSep = cellsize * 0.1;
    keepIdx = [true; ds > minSep];
    
    nRemoved = sum(~keepIdx);
    if nRemoved > 0
        fprintf('    Removed %d duplicate points (%.1f%% of original)\n', ...
            nRemoved, 100*nRemoved/n);
    end
    
    xo = xo(keepIdx);
    yo = yo(keepIdx);
end

end
