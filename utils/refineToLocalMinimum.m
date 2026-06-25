function [xr, yr, zr] = refineToLocalMinimum(x, y, DEM, radius)
%REFINETOLOCAMINIMUM Snap points to local MINIMUM (for trench/valley)
%
%   [xr, yr, zr] = refineToLocalMinimum(x, y, DEM, radius)
%
%   Part of S2S ModuleA - Source-to-Sink Geomorphometric Analysis
%
%   For each input point, searches within 'radius' pixels in the DEM
%   for the minimum elevation and snaps the point there.
%
%   INPUTS:
%       x, y   - Coordinate vectors (from coarse analysis)
%       DEM    - GRIDobj to sample from (typically full-resolution)
%       radius - Search radius in pixels [default: 3]
%
%   OUTPUTS:
%       xr, yr - Refined coordinates (snapped to local minimum)
%       zr     - Elevation at refined locations
%
%   REQUIRES:
%       TopoToolbox 2.x
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   O(1) coordinate conversion + bounds clamping

if nargin < 4
    radius = 3;
end

[xr, yr, zr] = refineToLocalExtremum(x, y, DEM, radius, 'min');

end


%% ========================================================================
% LOCAL HELPER FUNCTION
% =========================================================================

function [xr, yr, zr] = refineToLocalExtremum(x, y, DEM, radius, mode)
%REFINETOLOCAEXTREMUM Internal function for min/max refinement
%
%   Optimized: Uses O(1) coordinate-to-index conversion
%   instead of O(gridSize) min(abs()) search per point
%   Includes bounds clamping for edge cases.

x = x(:);
y = y(:);
n = numel(x);

xr = x;
yr = y;
zr = zeros(n, 1);

cellsize = DEM.cellsize;
[X, Y] = getcoordinates(DEM);
[nRows, nCols] = size(DEM.Z);

% Pre-compute coordinate conversion parameters (O(1) per point instead of O(gridSize))
x0 = X(1);
y0 = Y(1);
dx = abs(X(2) - X(1));
dy = abs(Y(2) - Y(1));  % Note: Y may be decreasing (north to south)
y_decreasing = Y(1) > Y(end);

for i = 1:n
    % Direct coordinate-to-index conversion: O(1) instead of O(gridSize)
    col0 = round((x(i) - x0) / dx) + 1;
    if y_decreasing
        row0 = round((y0 - y(i)) / dy) + 1;
    else
        row0 = round((y(i) - y0) / dy) + 1;
    end
    
    % CRITICAL: Clamp to valid range (prevents index-out-of-bounds)
    col0 = max(1, min(col0, nCols));
    row0 = max(1, min(row0, nRows));
    
    % Define search window
    rowMin = max(1, row0 - radius);
    rowMax = min(nRows, row0 + radius);
    colMin = max(1, col0 - radius);
    colMax = min(nCols, col0 + radius);
    
    % Extract window
    Z_window = DEM.Z(rowMin:rowMax, colMin:colMax);
    
    % Find extremum
    if strcmp(mode, 'min')
        [extremeZ, idx] = min(Z_window(:), [], 'omitnan');
    else
        [extremeZ, idx] = max(Z_window(:), [], 'omitnan');
    end
    
    if ~isnan(extremeZ)
        [localRow, localCol] = ind2sub(size(Z_window), idx);
        globalRow = rowMin + localRow - 1;
        globalCol = colMin + localCol - 1;
        
        xr(i) = X(globalCol);
        yr(i) = Y(globalRow);
        zr(i) = extremeZ;
    else
        % Fallback to interpolation
        zr(i) = interp(DEM, x(i), y(i));
    end
end

end
