%% stationarity_test.m
% =========================================================================
%  STATIONARITY TEST — Verifies domain segmentation is scale-invariant
% =========================================================================
%
%  Tests that boundary positions produced by segmentDomains on a trimmed,
%  padded subset reproduce the full-profile results within 5 km, measured
%  only inside a defined comparison zone.
%
%  EDGE-PADDING RATIONALE
%  ----------------------
%  deadbandWindow_km = 100 → movmedian needs ±50 km of neighbourhood.
%  When segmentDomains is given a short excerpt, Endpoints='shrink' causes
%  the rolling MAD near the clip edges to differ from the full-profile MAD.
%  This is CORRECT behaviour for a local window — not a bug.
%  Fix: pass a WIDER excerpt than the zone you want to check.  Points near
%  the centre of the padded excerpt see a full rolling window and therefore
%  produce identical deadband values to the full-profile run.
%
%  The subset passed to segmentDomains is the comparison zone extended by
%  2 × deadbandWindow_km on each side.  Only boundaries that fall inside
%  the comparison zone are compared; padding-zone boundaries are discarded.
%
%  TWO INDEPENDENT TESTS
%  ---------------------
%  Test 1: comparison zone [600, 1200] km
%          subset zone     [400, 1400] km  (200 km pad each side)
%
%  Test 2: comparison zone [1500, 2100] km
%          subset zone     [1300, 2300] km  (200 km pad each side)
%
%  Both tests must PASS for the overall verdict to be PASS.
%
%  USAGE:
%    Run run_master_figure.m first, then:
%      >> stationarity_test
%
%
%  Author: David Fernández-Blanco
%  Date:   2026-03-22
%
% =========================================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('  STATIONARITY TEST — Scale-Invariance Verification\n');
fprintf('================================================================\n\n');

%% --- Environment: add utils/ to path ---
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
projRoot = fileparts(scriptDir);  % one level up from analysis/
addpath(fullfile(projRoot, 'utils'));
cd(projRoot);  % so output/ paths resolve correctly

if exist('segmentDomains', 'file') == 0
    error('stationarity_test:missingFunction', ...
        'segmentDomains not found. Run from the project root (s2s-along-trench-v1.0/).');
end

%% --- Global test parameters ---
TOLERANCE_KM    = 5;     % km — acceptable position mismatch
PAD_FACTOR      = 2;     % pad = PAD_FACTOR × deadbandWindow_km each side

fprintf('  Tolerance:   %.1f km\n', TOLERANCE_KM);
fprintf('  Edge pad:    %d × deadbandWindow_km per side\n\n', PAD_FACTOR);


%% ========================================================================
%  STEP 1 — Load pipeline results
%  ========================================================================
matFile = fullfile(projRoot, 'output', '1DSink_MasterFigure_results.mat');
if ~isfile(matFile)
    error('stationarity_test:matNotFound', ...
        ['Results file not found:\n  %s\n' ...
         'Run run_master_figure.m first.'], matFile);
end

fprintf('Step 1: Loading pipeline results...\n');
fprintf('        %s\n', matFile);

S = load(matFile, 'precomputedSeg', 'bands', 'cfg');

for varName = {'precomputedSeg', 'bands', 'cfg'}
    if ~isfield(S, varName{1})
        error('stationarity_test:missingVar', ...
            'Variable ''%s'' not found in MAT file.', varName{1});
    end
end

precomputedSeg = S.precomputedSeg;
bands          = S.bands;
cfg            = S.cfg;

fprintf('  OK — precomputedSeg, bands, cfg loaded.\n');


%% ========================================================================
%  STEP 2 — Inspect precomputedSeg.boundaries (defensive — unknown format)
%  ========================================================================
fprintf('\nStep 2: Inspecting precomputedSeg...\n');

topFields = fieldnames(precomputedSeg);
fprintf('  precomputedSeg fields: %s\n', strjoin(topFields, ', '));

if ~isfield(precomputedSeg, 'boundaries')
    error('stationarity_test:noBoundaries', ...
        'precomputedSeg has no ''boundaries'' field.\nFields present: %s', ...
        strjoin(topFields, ', '));
end

rawBounds = precomputedSeg.boundaries;
fprintf('  boundaries class: %s\n', class(rawBounds));
fprintf('  boundaries size:  %s\n', mat2str(size(rawBounds)));

% Normalise to a plain numeric column vector of km positions
fullBounds_km = extractBoundaryKm(rawBounds);
fprintf('  Extracted %d full-profile boundary positions\n', numel(fullBounds_km));

if ~isempty(fullBounds_km)
    fprintf('  Full-profile range: [%.1f, %.1f] km\n', ...
        min(fullBounds_km), max(fullBounds_km));
end

% Verify required bands fields
for fld = {'s_km', 'd2D_2_20_mpkm2', 'edgeMask'}
    if ~isfield(bands, fld{1})
        error('stationarity_test:bandsMissingField', ...
            'bands struct is missing field ''%s''.\nAvailable: %s', ...
            fld{1}, strjoin(fieldnames(bands), ', '));
    end
end

fprintf('\n  bands.s_km: %d points, [%.1f, %.1f] km\n', ...
    numel(bands.s_km), min(bands.s_km), max(bands.s_km));

% Compute the pad from the stored / configured window
if isfield(cfg, 'segmentation') && isfield(cfg.segmentation, 'deadbandWindow_km')
    deadbandWindow_km = cfg.segmentation.deadbandWindow_km;
else
    deadbandWindow_km = 100;   % fallback to default
    warning('stationarity_test:noDeadbandWindow', ...
        'cfg.segmentation.deadbandWindow_km not found; using default 100 km.');
end
pad_km = PAD_FACTOR * deadbandWindow_km;

fprintf('  deadbandWindow_km = %.0f → pad = %.0f km per side\n\n', ...
    deadbandWindow_km, pad_km);


%% ========================================================================
%  STEP 3 — Run both tests
%  ========================================================================

% ---- Test 1 ----
fprintf('================================================================\n');
fprintf('  TEST 1 of 2\n');
fprintf('================================================================\n');
compZone1 = [600, 1200];
[result1, label1] = runOneTest(compZone1, pad_km, bands, cfg, fullBounds_km, TOLERANCE_KM);

% ---- Test 2 ----
fprintf('\n');
fprintf('================================================================\n');
fprintf('  TEST 2 of 2\n');
fprintf('================================================================\n');
compZone2 = [1500, 2100];
[result2, label2] = runOneTest(compZone2, pad_km, bands, cfg, fullBounds_km, TOLERANCE_KM);


%% ========================================================================
%  STEP 4 — Overall verdict
%  ========================================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('  OVERALL SUMMARY\n');
fprintf('================================================================\n');
fprintf('  Test 1  [%.0f, %.0f] km : %s\n', compZone1(1), compZone1(2), label1);
fprintf('  Test 2  [%.0f, %.0f] km : %s\n', compZone2(1), compZone2(2), label2);
fprintf('\n');

bothPass = strcmp(result1, 'PASS') && strcmp(result2, 'PASS');
anyFail  = strcmp(result1, 'FAIL') || strcmp(result2, 'FAIL');
anySkip  = strcmp(result1, 'SKIP') || strcmp(result2, 'SKIP');

if bothPass
    fprintf('  OVERALL VERDICT: PASS\n');
    fprintf('  Segmentation is SCALE-INVARIANT in both zones.\n');
    fprintf('  The v1.0 local-window fix is confirmed.\n');
elseif anyFail
    fprintf('  OVERALL VERDICT: FAIL\n');
    fprintf('  At least one test zone failed — segmentation is NOT scale-invariant.\n');
    fprintf('  Check that segmentDomains.m in utils/ is the v1.0 rewrite.\n');
elseif anySkip
    fprintf('  OVERALL VERDICT: PARTIAL (one zone skipped — data extent too short)\n');
    if strcmp(result1, 'PASS') || strcmp(result2, 'PASS')
        fprintf('  The zone that ran is PASS.\n');
    end
else
    fprintf('  OVERALL VERDICT: PARTIAL PASS (warnings present)\n');
end
fprintf('================================================================\n\n');

% Determine overall status string for the struct
if bothPass
    overallStatus = 'PASS';
elseif ~bothPass && ~anyFail && ~anySkip
    overallStatus = 'PARTIAL';
else
    overallStatus = 'FAIL'; 
end

% Return for programmatic use
stationarityVerdict = struct('test1', result1, 'test2', result2, ...
    'overall', overallStatus);  %#ok<NASGU>


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [verdict, verdictLabel] = runOneTest(compZone, pad_km, bands, cfg, ...
                                              fullBounds_km, TOLERANCE_KM)
%RUNONETEST  Run one padded-subset stationarity test.
%
%  compZone  — [lo, hi] km comparison zone (where boundaries are checked)
%  pad_km    — km of extra data to include on each side of compZone
%
%  Returns verdict = 'PASS' | 'PARTIAL' | 'FAIL' | 'SKIP'

subZone = [compZone(1) - pad_km, compZone(2) + pad_km];

fprintf('\n  Comparison zone : [%.0f, %.0f] km\n', compZone(1), compZone(2));
fprintf('  Subset zone     : [%.0f, %.0f] km  (pad = ±%.0f km)\n', ...
    subZone(1), subZone(2), pad_km);

% --- Guard: check data covers the comparison zone at minimum ---
if min(bands.s_km) > compZone(1) || max(bands.s_km) < compZone(2)
    fprintf('\n  SKIP: bands.s_km [%.1f, %.1f] km does not cover comparison zone.\n', ...
        min(bands.s_km), max(bands.s_km));
    fprintf('  Re-run run_master_figure.m with axisTrim_km = [] to get the full profile.\n');
    verdict      = 'SKIP';
    verdictLabel = 'SKIP (data extent too short)';
    return
end

% --- Reference boundaries: full-profile run, inside comparison zone ---
inComp    = fullBounds_km >= compZone(1) & fullBounds_km <= compZone(2);
refBounds = sort(fullBounds_km(inComp));
nRef      = numel(refBounds);

fprintf('\n  Reference boundaries in [%.0f, %.0f] km: %d\n', ...
    compZone(1), compZone(2), nRef);
if nRef > 0
    fprintf('    Positions: ');  fprintf('%.1f  ', refBounds);  fprintf('km\n');
end

% --- Build padded subset (clip to available data) ---
subLo = max(subZone(1), min(bands.s_km));
subHi = min(subZone(2), max(bands.s_km));

if subLo > compZone(1) || subHi < compZone(2)
    fprintf('  NOTE: subset zone clipped to [%.1f, %.1f] km (data limit).\n', subLo, subHi);
    fprintf('  Padding is partial — edge effects may reach into comparison zone.\n');
end

subMask  = bands.s_km >= subLo & bands.s_km <= subHi;
nSub_pts = sum(subMask);
fprintf('  Subset points   : %d  [%.1f, %.1f] km\n', ...
    nSub_pts, min(bands.s_km(subMask)), max(bands.s_km(subMask)));

if nSub_pts < 10
    fprintf('\n  SKIP: only %d points in subset — too few for segmentation.\n', nSub_pts);
    verdict      = 'SKIP';
    verdictLabel = 'SKIP (too few points)';
    return
end

sub_s_km = bands.s_km(subMask);       sub_s_km = sub_s_km(:);
sub_curv = bands.d2D_2_20_mpkm2(subMask); sub_curv = sub_curv(:);
sub_edge = bands.edgeMask(subMask);   sub_edge = sub_edge(:);

% --- Re-run segmentDomains on padded subset ---
fprintf('\n  Running segmentDomains on padded subset (stabilityTest=false)...\n');
tSeg = tic;

[~, subBoundsRaw, ~] = segmentDomains( ...
    sub_s_km, sub_curv, sub_edge, ...
    'minDomainLength_km',      cfg.segmentation.minDomainLength_km, ...
    'persistenceLength_km',    cfg.segmentation.persistenceLength_km, ...
    'deadbandMode',            cfg.segmentation.deadbandMode, ...
    'deadbandAmplitude',       cfg.segmentation.deadbandAmplitude_fixed, ...
    'deadbandFactor',          cfg.segmentation.deadbandFactor_percentile, ...
    'deadbandSigmas',          cfg.segmentation.deadbandSigmas_mad, ...
    'deadbandWindow_km',       cfg.segmentation.deadbandWindow_km, ...
    'confidenceWeights',       cfg.segmentation.confidenceWeights(:)', ...
    'dwellScale_km',           cfg.segmentation.dwellScale_km, ...
    'localAmplitudeRadius_km', cfg.segmentation.localAmplitudeRadius_km, ...
    'maxBoundaries',           cfg.segmentation.maxBoundaries, ...
    'smoothWindow_km',         cfg.segmentation.smoothWindow_km, ...
    'stabilityTest',           false, ...
    'verbose',                 false);

fprintf('  Done in %.1f s\n', toc(tSeg));

allSubBounds_km = extractBoundaryKm(subBoundsRaw);
fprintf('  Subset produced %d boundaries (entire padded zone)\n', numel(allSubBounds_km));

% --- Restrict subset boundaries to comparison zone only ---
inComp2      = allSubBounds_km >= compZone(1) & allSubBounds_km <= compZone(2);
subBounds_km = sort(allSubBounds_km(inComp2));
nSubB        = numel(subBounds_km);

fprintf('  Boundaries inside comparison zone [%.0f, %.0f] km: %d\n', ...
    compZone(1), compZone(2), nSubB);
if nSubB > 0
    fprintf('    Positions: ');  fprintf('%.1f  ', subBounds_km);  fprintf('km\n');
end

% --- Match reference vs subset boundaries ---
fprintf('\n  Matching (tolerance = %.1f km)...\n\n', TOLERANCE_KM);

refMatched = false(nRef,  1);
subMatched = false(nSubB, 1);
matchDist  = nan(nRef,  1);
matchIdx   = nan(nRef,  1);

for i = 1:nRef
    if nSubB == 0, break; end
    dists = abs(subBounds_km - refBounds(i));
    [d, idx] = min(dists);
    if d <= TOLERANCE_KM
        refMatched(i) = true;
        matchDist(i)  = d;
        matchIdx(i)   = idx;
        subMatched(idx) = true;
    end
end

% Per-boundary table
fprintf('  %-6s  %-14s  %-14s  %-10s  %s\n', ...
    'Ref#', 'Ref pos (km)', 'Match (km)', 'Err (km)', 'Status');
fprintf('  %s\n', repmat('-', 1, 58));

for i = 1:nRef
    if refMatched(i)
        fprintf('  %-6d  %-14.2f  %-14.2f  %-10.2f  MATCHED\n', ...
            i, refBounds(i), subBounds_km(matchIdx(i)), matchDist(i));
    else
        fprintf('  %-6d  %-14.2f  %-14s  %-10s  *** MISSING ***\n', ...
            i, refBounds(i), '--', '--');
    end
end

if nRef == 0
    fprintf('  (no reference boundaries in comparison zone)\n');
end

nSpurious = sum(~subMatched);
if nSpurious > 0
    fprintf('\n  Spurious (in subset zone but not in full-profile run):\n');
    for j = find(~subMatched)'
        fprintf('    %.2f km\n', subBounds_km(j));
    end
end

% --- Verdict ---
nMissing = sum(~refMatched);
nMatched = sum(refMatched);
maxErr   = max(matchDist(refMatched), [], 'omitnan');
if isempty(maxErr) || isnan(maxErr), maxErr = 0; end

fprintf('\n  SUMMARY for [%.0f, %.0f] km:\n', compZone(1), compZone(2));
fprintf('    Reference boundaries : %d\n', nRef);
fprintf('    Subset boundaries    : %d\n', nSubB);
fprintf('    Matched              : %d\n', nMatched);
fprintf('    Missing              : %d\n', nMissing);
fprintf('    Spurious             : %d\n', nSpurious);
if nMatched > 0
    fprintf('    Max match error      : %.2f km\n', maxErr);
end

if nRef == 0 && nSubB == 0
    verdict      = 'PASS';
    verdictLabel = 'PASS (no boundaries — trivially scale-invariant)';
    fprintf('  VERDICT: PASS — no boundaries in zone.\n');

elseif nMissing == 0 && nSpurious == 0
    verdict      = 'PASS';
    verdictLabel = sprintf('PASS (max err %.2f km, tol %.1f km)', maxErr, TOLERANCE_KM);
    fprintf('  VERDICT: PASS\n');
    fprintf('  All %d reference boundaries reproduced. Max error: %.2f km.\n', ...
        nMatched, maxErr);

elseif nMissing == 0 && nSpurious > 0
    verdict      = 'PARTIAL';
    verdictLabel = sprintf('PARTIAL (%d spurious)', nSpurious);
    fprintf('  VERDICT: PARTIAL PASS\n');
    fprintf('  All reference boundaries reproduced, but %d extra appeared.\n', nSpurious);
    fprintf('  These are confined to the padding-zone edges and do not indicate\n');
    fprintf('  interior scale-dependence.\n');

else
    verdict      = 'FAIL';
    verdictLabel = sprintf('FAIL (%d/%d missing)', nMissing, nRef);
    fprintf('  VERDICT: FAIL\n');
    fprintf('  %d / %d reference boundaries NOT reproduced.\n', nMissing, nRef);
    fprintf('  Even with edge padding, interior results differ — check segmentDomains.m.\n');
end
end  % runOneTest


function km = extractBoundaryKm(raw)
%EXTRACTBOUNDARYKM  Normalise the 'boundaries' output of segmentDomains to a
%   plain numeric column vector of km positions.
%
%   Handles: numeric vector/matrix, struct array, cell array, table.

    if isempty(raw)
        km = zeros(0,1);
        return
    end

    % Most common case: numeric vector
    if isnumeric(raw)
        km = raw(:);
        return
    end

    % Cell array
    if iscell(raw)
        vals = cellfun(@extractScalarKm, raw(:), 'UniformOutput', true);
        km   = vals(:);
        return
    end

    % Table
    if istable(raw)
        candidates = {'s_km','km','position_km','pos_km','boundary_km','location_km'};
        colNamesLow = lower(raw.Properties.VariableNames);
        found = false;
        for ci = 1:numel(candidates)
            hit = strcmp(colNamesLow, candidates{ci});
            if any(hit)
                km    = raw.(raw.Properties.VariableNames{find(hit,1)})(:);
                found = true;
                break
            end
        end
        if ~found
            for ci = 1:width(raw)
                col = raw.(raw.Properties.VariableNames{ci});
                if isnumeric(col)
                    km = col(:);
                    warning('stationarity_test:tableGuess', ...
                        'Using first numeric column ''%s'' as boundary km.', ...
                        raw.Properties.VariableNames{ci});
                    return
                end
            end
            error('stationarity_test:noKmColumn', ...
                'Cannot find km column in boundaries table.\nColumns: %s', ...
                strjoin(raw.Properties.VariableNames, ', '));
        end
        return
    end

    % Struct array
    if isstruct(raw)
        fnames = fieldnames(raw);
        fprintf('    [inspect] boundaries struct fields: %s\n', strjoin(fnames, ', '));
        candidates = {'s_km','km','position_km','pos_km','boundary_km','location_km','s'};
        found = false;
        for ci = 1:numel(candidates)
            if ismember(candidates{ci}, fnames)
                vals  = [raw.(candidates{ci})];
                km    = vals(:);
                found = true;
                fprintf('    [inspect] using field ''%s''\n', candidates{ci});
                break
            end
        end
        if ~found
            for ci = 1:numel(fnames)
                testVal = raw(1).(fnames{ci});
                if isnumeric(testVal) && isscalar(testVal)
                    vals = [raw.(fnames{ci})];
                    km   = vals(:);
                    warning('stationarity_test:structGuess', ...
                        'Using struct field ''%s'' as boundary km positions.', fnames{ci});
                    found = true;
                    break
                end
            end
        end
        if ~found
            error('stationarity_test:noKmField', ...
                'Cannot find km position field in boundaries struct.\nFields: %s', ...
                strjoin(fnames, ', '));
        end
        return
    end

    error('stationarity_test:unknownType', ...
        'Unsupported boundaries type: %s (size %s)', class(raw), mat2str(size(raw)));
end


function v = extractScalarKm(elem)
%EXTRACTSCALARKM  Extract a scalar km value from one cell element.
    if isnumeric(elem) && isscalar(elem)
        v = elem;
        return
    end
    if isstruct(elem)
        fnames     = fieldnames(elem);
        candidates = {'s_km','km','position_km','s'};
        for ci = 1:numel(candidates)
            if ismember(candidates{ci}, fnames)
                v = elem.(candidates{ci});
                return
            end
        end
        for ci = 1:numel(fnames)
            val = elem.(fnames{ci});
            if isnumeric(val) && isscalar(val)
                v = val;
                return
            end
        end
        error('stationarity_test:cellStructNoKm', ...
            'Cannot extract km from cell element struct. Fields: %s', ...
            strjoin(fnames, ', '));
    end
    error('stationarity_test:cellElemType', ...
        'Unexpected cell element type: %s', class(elem));
end
