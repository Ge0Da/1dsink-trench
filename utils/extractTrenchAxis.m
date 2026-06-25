function axis = extractTrenchAxis(DEM, opts)
%EXTRACTTRENCHAXIS  Extract trench axis from DEM via drainage divide method
%
%   axis = extractTrenchAxis(DEM)
%   axis = extractTrenchAxis(DEM, 'coarseRes_m', 500, ...)
%
%   Extracts the trench axis from a bathymetric DEM by:
%     1. Coarsening the DEM for fast processing.
%     2. Inverting the DEM (trench becomes the highest ridge).
%     3. Computing flow routing on the inverted DEM.
%     4. Extracting the main drainage divide (= trench axis).
%     5. Ordering divide points into a continuous polyline.
%     6. Refining each point to the local depth minimum on the full-res DEM.
%     7. Smoothing the axis planform to remove grid staircase artefacts.
%     8. Optionally trimming to a specified along-axis range.
%
%   CONCEPT:
%     When a bathymetric DEM is inverted (multiplied by -1), the trench —
%     being the deepest feature — becomes the HIGHEST ridge. This ridge is
%     the main drainage divide in the inverted topography. Extracting it
%     gives the trench axis automatically, with no baseline assumptions.
%
%   INPUTS:
%     DEM   — GRIDobj (TopoToolbox) at full resolution, non-inverted.
%             Depth values should be negative (or at least: trench is the
%             minimum). NaN and land (z ≥ 0) are handled gracefully.
%
%   NAME-VALUE OPTIONS:
%     'coarseRes_m'           — Coarse-phase DEM resolution (m).       Default: 500
%     'minAreaKm2'            — Min drainage area for streams (km²).   Default: 100
%     'localRefineRadius_cells' — Snap radius at full res (pixels).   Default: 3
%     'smoothWindow_samples'  — Planform smoothing window (odd int).   Default: 11
%     'axisTrim_km'           — [start, end] along-axis range, or [].  Default: []
%     'refineToFullRes'       — Refine to full-res DEM (logical).      Default: true
%     'verbose'               — Print progress messages (logical).     Default: true
%
%   OUTPUTS:
%     axis.x       — X coordinates of axis (full-res units)
%     axis.y       — Y coordinates of axis (full-res units)
%     axis.z       — Depth at each axis point (from DEM)
%     axis.s_m     — Arc length along axis (metres), starting at 0
%     axis.s_km    — Arc length along axis (km)
%     axis.nNodes  — Number of axis points
%     axis.length_m  — Total axis length (m)
%     axis.length_km — Total axis length (km)
%     axis.depthRange — [min, max] depth along axis
%     axis.coarseRes  — Actual coarse resolution used (m)
%     axis.fullRes    — Native DEM resolution (m)
%     axis.params     — Structure of all parameters used
%
%   DEPENDENCIES:
%     TopoToolbox (GRIDobj, FLOWobj, STREAMobj)
%     src/axis/ utilities:
%       resampleToCoarse, computeFlowRouting, extractMainDivide,
%       orderPointsByConnectivity, refineToLocalMinimum,
%       computeArcLength, sampleFromDEM
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
    DEM
    opts.coarseRes_m            (1,1) double {mustBePositive} = 500
    opts.minAreaKm2             (1,1) double {mustBePositive} = 100
    opts.localRefineRadius_cells (1,1) double {mustBeNonnegative} = 3
    opts.smoothWindow_samples   (1,1) double {mustBePositive} = 11
    opts.axisTrim_km                                          = []
    opts.refineToFullRes        (1,1) logical                 = true
    opts.verbose                (1,1) logical                 = true
end

% Enforce odd smoothing window (even windows are asymmetric)
if mod(opts.smoothWindow_samples, 2) == 0
    opts.smoothWindow_samples = opts.smoothWindow_samples + 1;
    if opts.verbose
        fprintf('  [extractTrenchAxis] smoothWindow rounded up to %d (must be odd)\n', ...
            opts.smoothWindow_samples);
    end
end

% Validate axisTrim_km format
if ~isempty(opts.axisTrim_km)
    if ~isnumeric(opts.axisTrim_km) || numel(opts.axisTrim_km) ~= 2
        error('extractTrenchAxis:badTrim', ...
            'axisTrim_km must be [] or a 1x2 vector [start_km, end_km].');
    end
    if opts.axisTrim_km(1) >= opts.axisTrim_km(2)
        error('extractTrenchAxis:badTrim', ...
            'axisTrim_km(1) must be < axisTrim_km(2). Got [%.1f, %.1f].', ...
            opts.axisTrim_km(1), opts.axisTrim_km(2));
    end
end

% Store the full-resolution DEM separately (used for refinement later)
DEM_fullres = DEM;

if opts.verbose
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  TRENCH AXIS EXTRACTION (drainage divide on inverted DEM)\n');
    fprintf('================================================================\n');
    fprintf('  DEM resolution:         %.0f m\n', DEM.cellsize);
    fprintf('  Coarse resolution:      %.0f m\n', opts.coarseRes_m);
    fprintf('  Min drainage area:      %.0f km²\n', opts.minAreaKm2);
    fprintf('  Refine radius:          %d pixels\n', opts.localRefineRadius_cells);
    fprintf('  Smooth window:          %d points\n', opts.smoothWindow_samples);
    if ~isempty(opts.axisTrim_km)
        fprintf('  Axis trim:              [%.1f, %.1f] km\n', ...
            opts.axisTrim_km(1), opts.axisTrim_km(2));
    end
    fprintf('================================================================\n\n');
end


%% ====================================================================
%  PHASE 1: COARSEN DEM
%  ====================================================================
if opts.verbose, fprintf('--- Phase 1: Resample to coarse resolution ---\n'); end

[DEM_coarse, resInfo] = resampleToCoarse(DEM, opts.coarseRes_m, opts.verbose);


%% ====================================================================
%  PHASE 2: INVERT DEM
%  ====================================================================
%  When multiplied by -1, the deepest trench becomes the highest ridge.
%  The drainage divide of this inverted surface IS the trench axis.

if opts.verbose, fprintf('\n--- Phase 2: Invert DEM ---\n'); end

DEM_inv = DEM_coarse;
DEM_inv.Z = -DEM_coarse.Z;

if opts.verbose
    fprintf('  Depth range (original): %.0f to %.0f m\n', ...
        min(DEM_coarse.Z(:), [], 'omitnan'), max(DEM_coarse.Z(:), [], 'omitnan'));
end


%% ====================================================================
%  PHASE 3: FLOW ROUTING
%  ====================================================================
if opts.verbose, fprintf('\n--- Phase 3: Flow routing on inverted DEM ---\n'); end

[FD, ~, ~] = computeFlowRouting(DEM_inv, opts.verbose);


%% ====================================================================
%  PHASE 4: STREAM NETWORK
%  ====================================================================
if opts.verbose, fprintf('\n--- Phase 4: Stream network ---\n'); end

% Convert km² to cells at the coarse resolution
minAreaCells = (opts.minAreaKm2 * 1e6) / (resInfo.coarseRes^2);

if opts.verbose
    fprintf('  Threshold: %.0f km² = %.0f cells at %.0f m\n', ...
        opts.minAreaKm2, minAreaCells, resInfo.coarseRes);
end

S = STREAMobj(FD, 'minarea', minAreaCells);

if opts.verbose
    fprintf('  Stream network: %d nodes\n', numel(S.x));
end


%% ====================================================================
%  PHASE 5: EXTRACT DRAINAGE DIVIDE
%  ====================================================================
if opts.verbose, fprintf('\n--- Phase 5: Extract main drainage divide ---\n'); end

[div_x, div_y] = extractMainDivide(FD, S, DEM_coarse, opts.verbose);


%% ====================================================================
%  PHASE 6: ORDER POINTS INTO POLYLINE
%  ====================================================================
if opts.verbose, fprintf('\n--- Phase 6: Order points along divide ---\n'); end

[axis_x, axis_y] = orderPointsByConnectivity(div_x, div_y, resInfo.coarseRes);

if opts.verbose
    fprintf('  Ordered axis: %d points\n', numel(axis_x));
end


%% ====================================================================
%  PHASE 7: REFINE TO FULL RESOLUTION
%  ====================================================================
if opts.verbose, fprintf('\n--- Phase 7: Refine to full-resolution DEM ---\n'); end

if opts.refineToFullRes && resInfo.scaleFactor > 1 && opts.localRefineRadius_cells > 0
    % Snap each coarse-phase point to the local depth minimum on the
    % full-resolution DEM within a radius of localRefineRadius_cells.
    if opts.verbose
        fprintf('  Snapping to local minimum (radius: %d pixels)...\n', ...
            opts.localRefineRadius_cells);
    end
    [axis_x, axis_y, axis_z] = refineToLocalMinimum( ...
        axis_x, axis_y, DEM_fullres, opts.localRefineRadius_cells);

    if opts.verbose
        fprintf('  Depths sampled from %.0f m resolution DEM\n', DEM_fullres.cellsize);
    end
else
    % No refinement — sample depth from whichever DEM is appropriate
    if resInfo.scaleFactor > 1 && opts.refineToFullRes
        axis_z = sampleFromDEM(axis_x, axis_y, DEM_fullres);
    else
        axis_z = sampleFromDEM(axis_x, axis_y, DEM_coarse);
    end
    if opts.verbose
        fprintf('  Refinement skipped (radius=0 or single resolution)\n');
    end
end


%% ====================================================================
%  PHASE 8: SMOOTH PLANFORM AND RE-SAMPLE DEPTH
%  ====================================================================
if opts.smoothWindow_samples > 1 && numel(axis_x) > opts.smoothWindow_samples
    if opts.verbose, fprintf('\n--- Phase 8: Smooth axis planform ---\n'); end

    % Smooth (x, y) to remove grid staircase artefacts
    axis_x = smoothdata(axis_x, 'movmean', opts.smoothWindow_samples);
    axis_y = smoothdata(axis_y, 'movmean', opts.smoothWindow_samples);

    % Re-sample depth from the DEM at the smoothed locations
    % (ensures z values are always true DEM values, not interpolated artefacts)
    axis_z = sampleFromDEM(axis_x, axis_y, DEM_fullres);

    % Also smooth z to remove pixel-scale noise before differentiation
    % (movmedian is robust to outliers: spikes, NaN artifacts)
    axis_z = smoothdata(axis_z, 'movmedian', opts.smoothWindow_samples);

    if opts.verbose
        fprintf('  Smoothed (x,y,z) with %d-point window\n', opts.smoothWindow_samples);
    end
end


%% ====================================================================
%  COMPUTE ARC LENGTH
%  ====================================================================
axis_s = computeArcLength(axis_x, axis_y);


%% ====================================================================
%  OPTIONAL TRIMMING
%  ====================================================================
if ~isempty(opts.axisTrim_km)
    if opts.verbose, fprintf('\n--- Trimming axis to [%.1f, %.1f] km ---\n', ...
            opts.axisTrim_km(1), opts.axisTrim_km(2)); end

    fullLength_km = axis_s(end) / 1000;
    s_km = axis_s / 1000;

    % Validate requested range against actual length
    if opts.axisTrim_km(1) > fullLength_km
        warning('extractTrenchAxis:trimOutOfRange', ...
            'Trim start (%.1f km) exceeds axis length (%.1f km). Using full axis.', ...
            opts.axisTrim_km(1), fullLength_km);
    elseif opts.axisTrim_km(2) > fullLength_km
        warning('extractTrenchAxis:trimOutOfRange', ...
            'Trim end (%.1f km) exceeds axis length (%.1f km). Clamping to %.1f km.', ...
            opts.axisTrim_km(2), fullLength_km, fullLength_km);
        opts.axisTrim_km(2) = fullLength_km;
    end

    % Find points within the requested range
    keepMask = (s_km >= opts.axisTrim_km(1)) & (s_km <= opts.axisTrim_km(2));

    if sum(keepMask) < 100
        warning('extractTrenchAxis:tooFewPoints', ...
            'Trimming leaves only %d points (need ≥ 100). Using full axis.', sum(keepMask));
    else
        axis_x = axis_x(keepMask);
        axis_y = axis_y(keepMask);
        axis_z = axis_z(keepMask);

        % Recompute arc length from the new starting point
        axis_s = computeArcLength(axis_x, axis_y);

        if opts.verbose
            fprintf('  Trimmed: %.1f km → %.1f km (%d points)\n', ...
                fullLength_km, axis_s(end)/1000, numel(axis_x));
        end
    end
end


%% ====================================================================
%  BUILD OUTPUT STRUCTURE
%  ====================================================================

% Ensure column vectors
axis_x = axis_x(:);
axis_y = axis_y(:);
axis_z = axis_z(:);
axis_s = axis_s(:);

axis = struct();
axis.x         = axis_x;
axis.y         = axis_y;
axis.z         = axis_z;
axis.s_m       = axis_s;
axis.s_km      = axis_s / 1000;
axis.nNodes    = numel(axis_x);
axis.length_m  = axis_s(end);
axis.length_km = axis_s(end) / 1000;
axis.depthRange = [min(axis_z, [], 'omitnan'), max(axis_z, [], 'omitnan')];
axis.coarseRes = resInfo.coarseRes;
axis.fullRes   = DEM_fullres.cellsize;
axis.method    = 'drainage-divide-inverted';

% Record parameters for reproducibility
axis.params = opts;
axis.params.effectiveMinAreaCells = minAreaCells;

if opts.verbose
    fprintf('\n================================================================\n');
    fprintf('  AXIS EXTRACTION COMPLETE\n');
    fprintf('================================================================\n');
    fprintf('  Nodes:        %d\n', axis.nNodes);
    fprintf('  Length:        %.1f km\n', axis.length_km);
    fprintf('  Depth range:   %.0f to %.0f m\n', axis.depthRange(1), axis.depthRange(2));
    fprintf('  Resolution:    %.0f m (coarse) → %.0f m (full)\n', ...
        axis.coarseRes, axis.fullRes);
    fprintf('================================================================\n\n');
end

end
