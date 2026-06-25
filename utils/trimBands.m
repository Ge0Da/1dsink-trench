function bandsTrim = trimBands(bands, trimRange_km)
%TRIMBANDS Clip a bands structure to a sub-range of along-axis distance
%
%   bandsTrim = trimBands(bands, [s_start_km, s_end_km])
%
%   PURPOSE:
%   ========
%   Trims all parallel arrays in a 'bands' structure (from computeTrenchBands)
%   to a contiguous sub-range of the along-trench axis.  This is used when
%   the user wants to DISPLAY a portion of the trench while running the
%   ANALYSIS (spectral decomposition, segmentation) on the full profile.
%
%   Trimming at the display stage — rather than before analysis — ensures
%   that window-dependent algorithms (segmentation, spectral leakage) are
%   not affected by the choice of display range.
%
%   INPUTS:
%       bands         - Structure from computeTrenchBands().
%       trimRange_km  - [s_start, s_end] in km.  Empty [] = no-op.
%
%   OUTPUTS:
%       bandsTrim     - Trimmed copy of 'bands' with updated nPoints
%                       and length_km.
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

    % No-op if trim range is empty
    if isempty(trimRange_km)
        bandsTrim = bands;
        return;
    end

    validateattributes(trimRange_km, {'numeric'}, {'numel', 2});
    assert(trimRange_km(1) < trimRange_km(2), ...
        'trimBands:badRange', 'trimRange_km(1) must be < trimRange_km(2)');

    s = bands.s_km(:);
    n = numel(s);
    keep = s >= trimRange_km(1) & s <= trimRange_km(2);

    if sum(keep) < 2
        warning('trimBands:tooFew', ...
            'Trim range [%.1f, %.1f] retains fewer than 2 points.', ...
            trimRange_km(1), trimRange_km(2));
    end

    % --- Trim all top-level numeric arrays of length n ---
    bandsTrim = bands;  % copy scalars and non-array fields
    topFields = fieldnames(bands);

    for i = 1:numel(topFields)
        f = topFields{i};
        val = bands.(f);

        if isnumeric(val) && numel(val) == n
            bandsTrim.(f) = val(keep);
        elseif islogical(val) && numel(val) == n
            bandsTrim.(f) = val(keep);
        elseif isstruct(val)
            % Recurse one level (handles bands.lowpass, bands.highpass, etc.)
            subFields = fieldnames(val);
            for j = 1:numel(subFields)
                sf = subFields{j};
                sv = val.(sf);
                if (isnumeric(sv) || islogical(sv)) && numel(sv) == n
                    bandsTrim.(f).(sf) = sv(keep);
                end
            end
        end
    end

    % --- Update metadata ---
    bandsTrim.nPoints   = sum(keep);
    bandsTrim.length_km = max(bandsTrim.s_km) - min(bandsTrim.s_km);
    % NOTE: s_km is NOT shifted to zero.  Absolute along-trench coordinates
    % are preserved so that pre-computed domain boundaries (also in absolute
    % coordinates) align without adjustment.
end
