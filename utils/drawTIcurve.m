function drawTIcurve(ax, bands, metrics, opts, segmentBounds, ...
    hasSamples, sample_s_km, sample_names, ...
    xLim, colorTI, colorSample, colorBoundary, ...
    titleFS, labelFS, tickFS, annotFS, ds)
%DRAWTICURVE Draw the composite Trapping Index (TI) banner into a caller axes.
%
%   Extracted from plotTrenchBands_1DSink so the same content can be used as
%   profile 'e' of the master figure or as panel 'a' of the published
%   Figure 3 (TI curve over the four contaminant scatters).
%
%   INPUTS
%       ax                axes handle to draw into
%       bands             trench-band struct (with s_km field)
%       metrics           sink-metrics struct (with TI_full field)
%       opts              forwarded plot options (sampleMarkerSize, showSampleNames)
%       segmentBounds     1xN vector of along-axis km of domain boundaries
%       hasSamples        logical
%       sample_s_km       1xK along-axis km of sample stations
%       sample_names      1xK cell array of station identifiers
%       xLim              [xmin xmax] along-axis display range, km
%       colorTI, colorSample, colorBoundary    1x3 RGB triplets
%       titleFS, labelFS, tickFS, annotFS      font sizes (points)
%       ds                downsample stride for plotting

axes(ax);  hold(ax, 'on');

s_km   = bands.s_km;
TI_ds  = metrics.TI_full(1:ds:end);
s_ds   = s_km(1:ds:end);

% Filled TI band above the 0.3 baseline.
fill(ax, [s_ds; flipud(s_ds)], [TI_ds; 0.3*ones(size(TI_ds))], colorTI, ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax, s_ds, TI_ds, '-', 'Color', colorTI, 'LineWidth', 1.5);

% Reference lines: median 0.5 and elevated-trapping 0.55.
yline(ax, 0.5,  'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(ax, 0.55, '--', 'Color', [0.8 0 0], 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(ax, xLim(1) + diff(xLim)*0.01, 0.57, 'Elevated trapping', ...
    'FontSize', annotFS, 'Color', [0.8 0 0], 'FontWeight', 'bold');

ylim(ax, [0.3, 0.7]);

% Domain boundaries.
drawSegmentLines(segmentBounds, ax, colorBoundary);

% Sample markers + labels.
if hasSamples
    sample_TI = interp1(s_km, metrics.TI_full, sample_s_km, 'linear', NaN);
    scatter(ax, sample_s_km, sample_TI, opts.sampleMarkerSize, colorSample, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if opts.showSampleNames
        plotSampleLabels(ax, sample_s_km, sample_TI, sample_names, 0.03, ...
            'OneSided', false, 'MinSep', 70);
    end
end

xlabel(ax, 'Distance along trench (km)', 'FontSize', labelFS);
ylabel(ax, 'Trapping Index', 'FontSize', labelFS);
xlim(ax, xLim);
set(ax, 'FontSize', tickFS);
box(ax, 'on');
grid(ax, 'on');

% Mirror tick marks on the right y-axis (numbers + ticks, no label).
yyaxis(ax, 'right');
ylim(ax, [0.3, 0.7]);
set(ax, 'YColor', [0.15 0.15 0.15], 'FontSize', tickFS);
ylabel(ax, '');
yyaxis(ax, 'left');

hold(ax, 'off');

end
