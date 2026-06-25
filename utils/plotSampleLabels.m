function plotSampleLabels(ax, s_km, y_vals, names, offset, varargin)
%PLOTSAMPLELABELS Label samples with white-box text and smart collision avoidance.
%
%   Uses a single text object per sample with BackgroundColor for a clean,
%   perfectly centered white box behind each label.  Consecutive samples
%   closer than `MinSep` km are staggered vertically across up to six tiers
%   so the labels do not overlap.  Optional `'OneSided'` mode keeps every
%   label above the data point on screen, which is required for depth
%   panels (where YDir is reversed) and helpful for tightly clustered
%   datasets.
%
%   Syntax:
%       plotSampleLabels(ax, s_km, y_vals, names, offset)
%       plotSampleLabels(..., 'OneSided', true)
%       plotSampleLabels(..., 'FontSize', 8)
%       plotSampleLabels(..., 'MinSep', 100)
%       plotSampleLabels(..., 'FontColor', [0 0 0])
%       plotSampleLabels(..., 'BgColor',   [1 1 1])
%
%   Inputs:
%       ax       axes handle
%       s_km     along-axis positions of samples (km)
%       y_vals   y-coordinates at which the sample markers sit
%       names    cell array of label strings (one per sample)
%       offset   base vertical shift from the data point, in data units
%
%   Name-Value Options:
%       'OneSided'  (false)        true keeps every label above the data
%                                  point on screen
%       'FontSize'  (7)            label font size in points
%       'MinSep'    (80)           samples closer than this in km get
%                                  staggered onto adjacent tiers
%       'FontColor' ([0 0 0])      text colour
%       'BgColor'   ([1 1 1])      background box colour

    p = inputParser;
    addParameter(p, 'OneSided',  false);
    addParameter(p, 'FontSize',  7);
    addParameter(p, 'MinSep',    80);
    addParameter(p, 'FontColor', [0 0 0]);
    addParameter(p, 'BgColor',   [1 1 1]);
    parse(p, varargin{:});
    o = p.Results;

    s_km   = s_km(:);
    y_vals = y_vals(:);
    nSamp  = numel(s_km);
    if nSamp == 0, return; end

    isReversed = strcmp(get(ax, 'YDir'), 'reverse');

    if isReversed
        o.OneSided = true;
    end

    maxTiers = 6;
    if o.OneSided
        tierMult = [1, 2.3, 3.6, 4.9, 6.2, 7.5];
    else
        tierMult = [1, -1.3, 2.6, -2.6, 4.0, -4.0];
    end

    [~, sortIdx] = sort(s_km);
    tier = zeros(nSamp, 1);
    for i = 2:nSamp
        prevPos = s_km(sortIdx(i-1));
        currPos = s_km(sortIdx(i));
        if abs(currPos - prevPos) < o.MinSep
            tier(i) = mod(tier(i-1) + 1, maxTiers);
        else
            tier(i) = 0;
        end
    end
    tierOrig = zeros(nSamp, 1);
    tierOrig(sortIdx) = tier;

    yl     = ylim(ax);
    yRange = abs(yl(2) - yl(1));

    for i = 1:nSamp
        if isnan(y_vals(i)) || i > numel(names) || isempty(names{i})
            continue;
        end

        t    = tierOrig(i);
        mult = tierMult(mod(t, maxTiers) + 1);
        if isReversed
            y_label = y_vals(i) - mult * abs(offset);
        else
            y_label = y_vals(i) + mult * abs(offset);
        end

        margin  = yRange * 0.04;
        y_label = max(y_label, yl(1) + margin);
        y_label = min(y_label, yl(2) - margin);

        if (y_label >= y_vals(i) && ~isReversed) || ...
           (y_label <= y_vals(i) &&  isReversed)
            vAlign = 'bottom';
        else
            vAlign = 'top';
        end

        labelDist = abs(y_label - y_vals(i));
        if labelDist > abs(offset) * 1.5
            plot(ax, [s_km(i) s_km(i)], [y_vals(i) y_label], '-', ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 0.4, ...
                'HandleVisibility', 'off');
        end

        text(ax, s_km(i), y_label, names{i}, ...
            'FontSize',            o.FontSize, ...
            'FontWeight',          'bold', ...
            'Color',               o.FontColor, ...
            'BackgroundColor',     o.BgColor, ...
            'Margin',              0.5, ...
            'EdgeColor',           'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   vAlign, ...
            'Clipping',            'on', ...
            'HandleVisibility',    'off');
    end
end
