function [div_x, div_y] = extractMainDivide(FD, S, DEM, verbose)
%EXTRACTMAINDIVIDE Extract the main drainage divide (trench axis)
%
%   [div_x, div_y] = extractMainDivide(FD, S, DEM, verbose)
%
%   Part of S2S ModuleA - Source-to-Sink Geomorphometric Analysis
%
%   VERSION 8.1 - Robust outlet detection using max accumulation
%   
%   THEORY:
%   =======
%   The main divide separates basins flowing in OPPOSITE directions (~180° apart).
%   
%   For each basin, the outlet is the cell with MAXIMUM flow accumulation.
%   Flow direction = angle from basin centroid to outlet location.
%   
%   Basins are grouped by flow direction (~180° apart = opposite sides of divide).
%
%   INPUTS:
%       FD      - FLOWobj
%       S       - STREAMobj (unused, kept for API compatibility)
%       DEM     - GRIDobj (for coordinate conversion)
%       verbose - Print progress [default: true]
%
%   OUTPUTS:
%       div_x, div_y - Coordinates of main divide (unordered, LINEAR structure)
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Robust outlet detection

if nargin < 4
    verbose = true;
end

[nRows, nCols] = size(DEM.Z);

% =========================================================================
% STEP 1: Compute all natural drainage basins
% =========================================================================

if verbose
    fprintf('  Computing all drainage basins...\n');
end

DB = drainagebasins(FD);
basinIDs = unique(DB.Z(DB.Z > 0));
nBasins = numel(basinIDs);

if verbose
    fprintf('  Found %d natural drainage basins\n', nBasins);
end

if nBasins < 2
    warning('Only %d basin found. Cannot extract divide.', nBasins);
    div_x = [];
    div_y = [];
    return;
end

% =========================================================================
% STEP 2: For each basin, compute flow direction (centroid → max accumulation)
% =========================================================================

if verbose
    fprintf('  Computing flow direction for each basin...\n');
end

A = flowacc(FD);

basin_angle = zeros(nBasins, 1);  % Flow direction (radians)
basin_area = zeros(nBasins, 1);   % Area for weighting
valid_count = 0;

for i = 1:nBasins
    bid = basinIDs(i);
    basin_mask = (DB.Z == bid);
    basin_area(i) = sum(basin_mask(:));
    
    % Skip very small basins (less than 10 cells)
    if basin_area(i) < 10
        basin_angle(i) = NaN;
        continue;
    end
    
    % Basin centroid
    [rows_in_basin, cols_in_basin] = find(basin_mask);
    centroid_row = mean(rows_in_basin);
    centroid_col = mean(cols_in_basin);
    
    % Find outlet: cell with MAXIMUM flow accumulation in this basin
    % This is where all the water in the basin converges
    acc_in_basin = A.Z;
    acc_in_basin(~basin_mask) = -Inf;
    [max_acc, outlet_idx] = max(acc_in_basin(:));
    [outlet_row, outlet_col] = ind2sub([nRows, nCols], outlet_idx);
    
    % Flow direction: angle from centroid to outlet
    % This is the direction water flows in this basin
    dy = outlet_row - centroid_row;
    dx = outlet_col - centroid_col;
    
    % Skip if centroid equals outlet (degenerate case)
    if abs(dx) < 0.5 && abs(dy) < 0.5
        basin_angle(i) = NaN;
        continue;
    end
    
    basin_angle(i) = atan2(dy, dx);
    valid_count = valid_count + 1;
end

if verbose
    fprintf('  Valid basins with computable flow direction: %d\n', valid_count);
end

if valid_count < 2
    warning('Fewer than 2 valid basins. Cannot extract divide.');
    div_x = [];
    div_y = [];
    return;
end

% Remove basins with no valid angle
valid_mask = ~isnan(basin_angle);
basin_angle_valid = basin_angle(valid_mask);
basin_area_valid = basin_area(valid_mask);
basinIDs_valid = basinIDs(valid_mask);
nValid = sum(valid_mask);

% =========================================================================
% STEP 3: Find dominant flow axis using area-weighted circular statistics
% =========================================================================

if verbose
    fprintf('  Finding dominant flow axis...\n');
end

% Double the angles to convert directions to axes
doubled_angles = 2 * basin_angle_valid;

% Area-weighted circular mean of doubled angles
weights = basin_area_valid / sum(basin_area_valid);
mean_cos = sum(weights .* cos(doubled_angles));
mean_sin = sum(weights .* sin(doubled_angles));

% The dominant axis (halve the result to get back to original angle space)
dominant_axis = atan2(mean_sin, mean_cos) / 2;

if verbose
    fprintf('  Dominant flow axis: %.1f° (and %.1f°)\n', ...
        rad2deg(dominant_axis), rad2deg(mod(dominant_axis + 180, 360) - 180 + 180));
end

% =========================================================================
% STEP 4: Split basins into two groups based on flow direction
% =========================================================================

if verbose
    fprintf('  Splitting basins by flow direction...\n');
end

basin_group = zeros(nBasins, 1);  % 0 = invalid, 1 or 2 = group

for i = 1:nBasins
    if isnan(basin_angle(i))
        basin_group(i) = 0;
        continue;
    end
    
    % Angular difference from dominant axis
    delta = basin_angle(i) - dominant_axis;
    
    % Normalize to [-π, π]
    while delta > pi, delta = delta - 2*pi; end
    while delta < -pi, delta = delta + 2*pi; end
    
    % Group 1: within 90° of dominant axis
    % Group 2: within 90° of opposite direction
    if abs(delta) <= pi/2
        basin_group(i) = 1;
    else
        basin_group(i) = 2;
    end
end

% Compute statistics for each group
group1_mask = basin_group == 1;
group2_mask = basin_group == 2;
group1_area = sum(basin_area(group1_mask));
group2_area = sum(basin_area(group2_mask));
group1_count = sum(group1_mask);
group2_count = sum(group2_mask);

% Compute mean angle for each group
if group1_count > 0
    g1_angles = basin_angle(group1_mask);
    g1_weights = basin_area(group1_mask) / sum(basin_area(group1_mask));
    group1_mean = atan2(sum(g1_weights .* sin(g1_angles)), ...
                        sum(g1_weights .* cos(g1_angles)));
else
    group1_mean = 0;
end

if group2_count > 0
    g2_angles = basin_angle(group2_mask);
    g2_weights = basin_area(group2_mask) / sum(basin_area(group2_mask));
    group2_mean = atan2(sum(g2_weights .* sin(g2_angles)), ...
                        sum(g2_weights .* cos(g2_angles)));
else
    group2_mean = 0;
end

% Angular separation between groups
separation = abs(group1_mean - group2_mean);
if separation > pi, separation = 2*pi - separation; end

total_area = sum(basin_area);
if verbose
    fprintf('  Group 1: %d basins (%.1f%% area), mean direction: %.1f°\n', ...
        group1_count, 100*group1_area/total_area, rad2deg(group1_mean));
    fprintf('  Group 2: %d basins (%.1f%% area), mean direction: %.1f°\n', ...
        group2_count, 100*group2_area/total_area, rad2deg(group2_mean));
    fprintf('  Angular separation: %.1f° (ideal: 180°)\n', rad2deg(separation));
end

if separation < pi/2
    warning('Angular separation only %.1f°. Main divide may not be well-defined.', ...
        rad2deg(separation));
end

if group1_count == 0 || group2_count == 0
    warning('One group is empty. Cannot extract main divide.');
    div_x = [];
    div_y = [];
    return;
end

% =========================================================================
% STEP 5: Relabel grid into two mega-basins
% =========================================================================

if verbose
    fprintf('  Merging basins into two mega-basins...\n');
end

mega_basin = zeros(nRows, nCols);

for i = 1:nBasins
    if basin_group(i) > 0
        bid = basinIDs(i);
        mega_basin(DB.Z == bid) = basin_group(i);
    end
end

side1_area = sum(mega_basin(:) == 1);
side2_area = sum(mega_basin(:) == 2);

if verbose
    fprintf('  Mega-basin 1: %.1f%% of domain\n', 100*side1_area/(nRows*nCols));
    fprintf('  Mega-basin 2: %.1f%% of domain\n', 100*side2_area/(nRows*nCols));
end

if side1_area == 0 || side2_area == 0
    warning('One mega-basin is empty. Cannot extract main divide.');
    div_x = [];
    div_y = [];
    return;
end

% =========================================================================
% STEP 6: Extract boundary between the two mega-basins
% =========================================================================

if verbose
    fprintf('  Extracting main divide (mega-basin boundary)...\n');
end

% Pad array to handle edges
MB_padded = padarray(mega_basin, [1 1], 0, 'both');

center = MB_padded(2:end-1, 2:end-1);
left   = MB_padded(2:end-1, 1:end-2);
right  = MB_padded(2:end-1, 3:end);
up     = MB_padded(1:end-2, 2:end-1);
down   = MB_padded(3:end,   2:end-1);

% Boundary: cell in group 1 adjacent to group 2, or vice versa
boundary_mask = ((center == 1) & ((left == 2) | (right == 2) | (up == 2) | (down == 2))) | ...
                ((center == 2) & ((left == 1) | (right == 1) | (up == 1) | (down == 1)));

nBoundaryPixels = sum(boundary_mask(:));

if verbose
    fprintf('  Boundary pixels: %d\n', nBoundaryPixels);
end

if nBoundaryPixels == 0
    warning('No boundary found between mega-basins.');
    div_x = [];
    div_y = [];
    return;
end

% Take largest connected component
CC = bwconncomp(boundary_mask, 8);
if CC.NumObjects > 1
    componentSizes = cellfun(@numel, CC.PixelIdxList);
    [maxSize, maxIdx] = max(componentSizes);
    
    if verbose
        fprintf('  Note: %d boundary components, using largest (%d pixels, %.1f%%)\n', ...
            CC.NumObjects, maxSize, 100*maxSize/nBoundaryPixels);
    end
    
    [rows, cols] = ind2sub([nRows, nCols], CC.PixelIdxList{maxIdx});
else
    [rows, cols] = find(boundary_mask);
end

% Convert to coordinates
[X, Y] = getcoordinates(DEM);
div_x = X(cols)';
div_y = Y(rows)';

if verbose
    fprintf('  Main divide: %d points\n', numel(div_x));
    
    % Sanity check
    expectedOrder = max(nRows, nCols) * 2;
    if numel(div_x) > expectedOrder * 5
        warning('Divide has %d points (expected ~%d). Check configuration.', ...
            numel(div_x), expectedOrder);
    elseif numel(div_x) < 100
        warning('Divide has only %d points. Check basin configuration.', numel(div_x));
    end
end

end
