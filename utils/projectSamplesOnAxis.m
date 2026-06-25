function projected = projectSamplesOnAxis(sampleFile, axis, opts)
%PROJECTSAMPLESONAXIS  Project sediment samples onto the trench axis
%
%   projected = projectSamplesOnAxis(sampleFile, axis)
%   projected = projectSamplesOnAxis(sampleFile, axis, 'colX', 'Easting', ...)
%
%   For each sample, finds the closest point on the trench axis polyline
%   using segment-by-segment orthogonal projection. Returns along-axis
%   distance (s_km) and signed lateral offset (d_km).
%
%   ALGORITHM:
%     For each sample point P and each axis segment [V1, V2]:
%       1. Compute the tangent t̂ and normal n̂ of the segment.
%       2. Project P onto the segment line: t = (P − V1) · t̂
%       3. Clamp t to [0, segment_length] (closest point on segment).
%       4. Compute distance from P to the clamped point.
%       5. Record the (s, d, distance) for the segment that minimises distance.
%     The signed lateral offset d = (P − V1) · n̂, where n̂ points 90°
%     counterclockwise from the tangent.
%
%   INPUTS:
%     sampleFile — Path to CSV/Excel file, OR a MATLAB table
%     axis       — Structure from extractTrenchAxis (needs .x, .y, .s_m)
%
%   NAME-VALUE OPTIONS:
%     'colX'    — Column name for X coordinate   (default: 'x')
%     'colY'    — Column name for Y coordinate   (default: 'y')
%     'colName' — Column name for sample ID       (default: 'name')
%     'showNames' — Include names in output       (default: true)
%     'verbose'   — Print progress messages       (default: true)
%
%   OUTPUTS:
%     projected.s_km      — Along-axis distance for each sample (km)
%     projected.s_axis    — Along-axis distance (m)
%     projected.d_km      — Signed lateral offset (km), + = right of axis
%     projected.x         — Original sample X coordinates
%     projected.y         — Original sample Y coordinates
%     projected.name      — Sample names (cell array of strings)
%     projected.nSamples  — Number of samples
%     projected.dist_km   — Euclidean distance to nearest axis point (km)
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

%% ====================================================================
%  PARSE INPUTS
%  ====================================================================
arguments
    sampleFile
    axis       struct
    opts.colX    (1,:) char = 'x'
    opts.colY    (1,:) char = 'y'
    opts.colName (1,:) char = 'name'
    opts.showNames (1,1) logical = true
    opts.verbose   (1,1) logical = true
end

%% ====================================================================
%  LOAD SAMPLE TABLE
%  ====================================================================
if ischar(sampleFile) || isstring(sampleFile)
    % Load from file path
    samplePath = char(sampleFile);
    if ~isfile(samplePath)
        error('projectSamplesOnAxis:fileNotFound', ...
            'Sample file not found: %s', samplePath);
    end
    try
        samples = readtable(samplePath);
    catch ME
        error('projectSamplesOnAxis:readFailed', ...
            'Failed to read sample file "%s": %s', samplePath, ME.message);
    end
    if opts.verbose
        fprintf('  Loaded %d samples from %s\n', height(samples), samplePath);
    end
elseif istable(sampleFile)
    samples = sampleFile;
elseif isstruct(sampleFile)
    % Convert struct to table
    samples = struct2table(sampleFile);
else
    error('projectSamplesOnAxis:badInput', ...
        'sampleFile must be a file path, table, or structure.');
end

%% ====================================================================
%  FIND COORDINATE COLUMNS (with robust fallback search)
%  ====================================================================
colNames = samples.Properties.VariableNames;

sp_x = findAndExtractColumn(samples, colNames, opts.colX, ...
    {'x', 'X', 'Easting', 'easting', 'Longitude', 'longitude', 'lon', 'Lon'});
sp_y = findAndExtractColumn(samples, colNames, opts.colY, ...
    {'y', 'Y', 'Northing', 'northing', 'Latitude', 'latitude', 'lat', 'Lat'});

if isempty(sp_x) || isempty(sp_y)
    error('projectSamplesOnAxis:missingCoords', ...
        'Could not find coordinate columns. Available columns: %s\nSearched for X: %s\nSearched for Y: %s', ...
        strjoin(colNames, ', '), opts.colX, opts.colY);
end

sp_x = double(sp_x(:));
sp_y = double(sp_y(:));

%% ====================================================================
%  FIND NAME COLUMN (with fallback search)
%  ====================================================================
sp_names = findAndExtractColumn(samples, colNames, opts.colName, ...
    {'name', 'Name', 'OriginalID', 'ID', 'SampleID', 'Station', 'station', 'Sample'});

sp_names = normaliseNames(sp_names, numel(sp_x));

%% ====================================================================
%  VALIDATE AXIS STRUCTURE
%  ====================================================================
if ~isfield(axis, 'x') || ~isfield(axis, 'y') || ~isfield(axis, 's_m')
    error('projectSamplesOnAxis:badAxis', ...
        'axis structure must contain fields: x, y, s_m');
end

ax_x = axis.x(:);
ax_y = axis.y(:);
ax_s = axis.s_m(:);
nAxisPts = numel(ax_x);

if nAxisPts < 2
    error('projectSamplesOnAxis:tooFewAxisPoints', ...
        'Axis has only %d points. Need at least 2.', nAxisPts);
end

%% ====================================================================
%  ORTHOGONAL PROJECTION — segment by segment
%  ====================================================================
nSamples = numel(sp_x);
proj_s    = nan(nSamples, 1);   % Along-axis distance (m)
proj_d    = nan(nSamples, 1);   % Signed lateral offset (m)
proj_dist = nan(nSamples, 1);   % Euclidean distance to axis (m)

if opts.verbose
    fprintf('  Projecting %d samples onto %d-point axis...\n', nSamples, nAxisPts);
end

% Pre-compute segment properties for speed
seg_dx = diff(ax_x);
seg_dy = diff(ax_y);
seg_len = sqrt(seg_dx.^2 + seg_dy.^2);

% Unit tangents and normals (90° counterclockwise)
seg_tx = seg_dx ./ max(seg_len, 1e-12);
seg_ty = seg_dy ./ max(seg_len, 1e-12);
seg_nx = -seg_ty;
seg_ny =  seg_tx;

nSegs = nAxisPts - 1;

for i = 1:nSamples
    px = sp_x(i);
    py = sp_y(i);

    if isnan(px) || isnan(py)
        continue;
    end

    % Vector from each segment start to the sample point
    vx = px - ax_x(1:nSegs);
    vy = py - ax_y(1:nSegs);

    % Project onto tangent: how far along each segment
    t_proj = vx .* seg_tx + vy .* seg_ty;

    % Project onto normal: signed lateral offset
    n_proj = vx .* seg_nx + vy .* seg_ny;

    % Clamp t to segment bounds [0, seg_len]
    t_clamped = max(0, min(seg_len, t_proj));

    % Closest point on each segment
    cx = ax_x(1:nSegs) + t_clamped .* seg_tx;
    cy = ax_y(1:nSegs) + t_clamped .* seg_ty;

    % Distance from sample to closest point on each segment
    dists = sqrt((px - cx).^2 + (py - cy).^2);

    % Find the segment with minimum distance
    [minDist, bestSeg] = min(dists);

    % Record the projection results
    proj_s(i)    = ax_s(bestSeg) + t_clamped(bestSeg);
    proj_d(i)    = n_proj(bestSeg);
    proj_dist(i) = minDist;
end

%% ====================================================================
%  BUILD OUTPUT STRUCTURE
%  ====================================================================
projected = struct();
projected.x         = sp_x;
projected.y         = sp_y;
projected.name      = sp_names;
projected.s_km      = proj_s / 1000;
projected.s_axis    = proj_s;           % metres
projected.d_km      = proj_d / 1000;
projected.dist_km   = proj_dist / 1000;
projected.nSamples  = nSamples;

if opts.verbose
    fprintf('  Projection complete:\n');
    validMask = ~isnan(proj_s);
    fprintf('    Projected:   %d / %d samples\n', sum(validMask), nSamples);
    if any(validMask)
        fprintf('    Along-axis:  %.1f to %.1f km\n', ...
            min(proj_s(validMask))/1000, max(proj_s(validMask))/1000);
        fprintf('    Max lateral:  %.1f km\n', max(abs(proj_d(validMask)))/1000);
    end
end

end


%% ====================================================================
%  LOCAL HELPER FUNCTIONS
%  ====================================================================

function val = findAndExtractColumn(tbl, colNames, primary, fallbacks)
%FINDANDEXTRACTCOLUMN  Search for a column by primary name, then fallbacks
%
%   Tries exact match on 'primary' first, then tries each fallback in
%   order. Returns the column data if found, or [] if not found.

    % Try primary name first
    if ismember(primary, colNames)
        val = tbl.(primary);
        return;
    end

    % Try each fallback in order
    for f = 1:numel(fallbacks)
        if ismember(fallbacks{f}, colNames)
            val = tbl.(fallbacks{f});
            return;
        end
    end

    % Not found
    val = [];
end


function nameCell = normaliseNames(names, nExpected)
%NORMALISENAMES  Convert various name formats to a cell array of strings
%
%   Handles: string, cell, char, numeric, or empty.
%   Falls back to 'S1', 'S2', ... if nothing is found.

    if isempty(names)
        nameCell = arrayfun(@(i) sprintf('S%d', i), (1:nExpected)', 'UniformOutput', false);
        return;
    end

    if isstring(names)
        nameCell = cellstr(names);
    elseif iscell(names)
        nameCell = names;
    elseif ischar(names)
        nameCell = {names};
    elseif isnumeric(names)
        nameCell = arrayfun(@(x) num2str(x), names(:), 'UniformOutput', false);
    elseif iscategorical(names)
        nameCell = cellstr(names);
    else
        nameCell = arrayfun(@(i) sprintf('S%d', i), (1:nExpected)', 'UniformOutput', false);
    end

    % Ensure column vector
    nameCell = nameCell(:);

    % Pad or trim to expected length
    if numel(nameCell) < nExpected
        extra = arrayfun(@(i) sprintf('S%d', i), ...
            (numel(nameCell)+1 : nExpected)', 'UniformOutput', false);
        nameCell = [nameCell; extra];
    elseif numel(nameCell) > nExpected
        nameCell = nameCell(1:nExpected);
    end
end
