function savedFiles = saveFigure(fig, outputDir, filename, formats, opts)
%SAVEFIGURE Save figure with consistent formatting and multiple output formats
%
%   saveFigure(fig, outputDir, filename)
%   saveFigure(fig, outputDir, filename, formats)
%   saveFigure(fig, outputDir, filename, formats, 'resolution', 300)
%   savedFiles = saveFigure(...)
%
%   INPUTS:
%       fig        - Figure handle or [] for current figure
%       outputDir  - Output directory path
%       filename   - Filename without extension
%       formats    - Single format string or cell array of formats
%                    Supported: 'png', 'pdf', 'svg', 'eps', 'fig', 'tif'
%                    Default: {'png', 'svg', 'pdf'}
%
%   NAME-VALUE OPTIONS:
%       'resolution'    - DPI for raster formats (default: 300)
%       'transparent'   - true/false for transparent background (default: true)
%       'verbose'       - Print status messages (default: true)
%
%   OUTPUTS:
%       savedFiles - Cell array of saved file paths
%
%   TRANSPARENCY NOTES:
%       PNG transparency requires both figure and axes Color set to 'none'.
%       PDF/SVG do not reliably support transparency; white background used.
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 3.0

arguments
    fig = []
    outputDir (1,:) char = '.'
    filename (1,:) char = 'figure'
    formats = {'png', 'svg', 'pdf'}
    opts.resolution (1,1) double = 600
    opts.transparent (1,1) logical = true
    opts.verbose (1,1) logical = true
end

%% Get figure handle
if isempty(fig)
    fig = gcf;
end

%% Ensure formats is cell array
if ischar(formats) || isstring(formats)
    formats = {char(formats)};
end

%% Ensure output directory exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Store original figure properties for restoration
origColor = get(fig, 'Color');
origInvertHardcopy = get(fig, 'InvertHardcopy');

% Store original axes colors
allAxes = findall(fig, 'Type', 'axes');
origAxesColors = cell(length(allAxes), 1);
for i = 1:length(allAxes)
    origAxesColors{i} = get(allAxes(i), 'Color');
end

%% Suppress common export warnings
warnState = warning('query', 'all');
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
warning('off', 'MATLAB:exportgraphics:CannotUseTransparency');
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');

%% Save in each format
savedFiles = {};

for i = 1:length(formats)
    fmt = lower(formats{i});
    fullPath = fullfile(outputDir, [filename, '.', fmt]);
    
    try
        switch fmt
            case 'png'
                if opts.transparent
                    % Set figure and axes to transparent BEFORE export
                    set(fig, 'Color', 'none');
                    set(fig, 'InvertHardcopy', 'off');
                    for j = 1:length(allAxes)
                        if isequal(origAxesColors{j}, [1 1 1])
                            set(allAxes(j), 'Color', 'none');
                        end
                    end
                    % Export with transparency
                    exportgraphics(fig, fullPath, 'Resolution', opts.resolution, ...
                        'BackgroundColor', 'none');
                    % Restore immediately
                    if isvalid(fig)
                        set(fig, 'Color', origColor);
                    end
                    set(fig, 'InvertHardcopy', origInvertHardcopy);
                    for j = 1:length(allAxes)
                        set(allAxes(j), 'Color', origAxesColors{j});
                    end
                else
                    print(fig, fullPath, '-dpng', sprintf('-r%d', opts.resolution));
                end
                
            case 'pdf'
                % Vector PDF (white background - transparency not reliable)
                exportgraphics(fig, fullPath, 'ContentType', 'vector', ...
                    'BackgroundColor', 'white');
                
            case 'svg'
                % SVG vector format
                % Note: print -dsvg can invalidate fig handle on some platforms.
                % Use exportgraphics where available (R2020a+), fall back to print.
                try
                    exportgraphics(fig, fullPath, 'ContentType', 'vector', ...
                        'BackgroundColor', 'white');
                catch
                    if isvalid(fig)
                        print(fig, fullPath, '-dsvg');
                    end
                end
                
            case 'eps'
                % EPS vector format
                print(fig, fullPath, '-depsc2');
                
            case 'tif'
                % TIFF raster
                print(fig, fullPath, '-dtiff', sprintf('-r%d', opts.resolution));
                
            case 'fig'
                % MATLAB figure file
                savefig(fig, fullPath);
                
            otherwise
                warning('saveFigure:unknownFormat', ...
                    'Unknown format: %s. Skipping.', fmt);
                continue;
        end
        
        savedFiles{end+1} = fullPath; %#ok<AGROW>
        
        if opts.verbose
            fprintf('  Saved: %s\n', fullPath);
        end
        
    catch ME
        warning('saveFigure:saveError', ...
            'Failed to save %s: %s', fullPath, ME.message);
    end
end

%% Restore warning state
warning(warnState);

%% Restore original figure properties (safety net)
if isvalid(fig)
    set(fig, 'Color', origColor);
    set(fig, 'InvertHardcopy', origInvertHardcopy);
    for i = 1:length(allAxes)
        if isvalid(allAxes(i))
            set(allAxes(i), 'Color', origAxesColors{i});
        end
    end
end

%% Summary
if opts.verbose && length(savedFiles) > 1
    fprintf('  → %d files saved for: %s\n', length(savedFiles), filename);
end

end
