function drawSegmentLines(bounds, ax, rgb)
%DRAWSEGMENTLINES Draw vertical dashed lines at segment boundaries.
%   Lines are drawn in a faint version of `rgb` and excluded from the legend.

    if isempty(bounds)
        return;
    end

    rgb_faint = rgb * 0.6 + 0.4;

    for b = bounds(:)'
        xline(ax, b, '--', 'Color', rgb_faint, 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
end
