function [DEM_coarse, info] = resampleToCoarse(DEM, targetRes, verbose)
%RESAMPLETOCOARSE Resample DEM to coarse resolution for fast processing
%
%   [DEM_coarse, info] = resampleToCoarse(DEM, targetRes, verbose)
%
%   Part of S2S ModuleA - Source-to-Sink Geomorphometric Analysis
%   Dual-Resolution Architecture: Coarse for speed, fine for accuracy
%
%   INPUTS:
%       DEM       - GRIDobj at native resolution
%       targetRes - Target coarse resolution in meters [default: 500]
%       verbose   - Print progress messages [default: true]
%
%   OUTPUTS:
%       DEM_coarse - Resampled GRIDobj (or original if already coarse)
%       info       - Structure with:
%                    .nativeRes    - Original resolution
%                    .coarseRes    - Actual coarse resolution
%                    .scaleFactor  - Resampling factor
%                    .speedup      - Approximate speedup (scaleFactor^2)
%                    .nativeCells  - Original cell count
%                    .coarseCells  - Coarse cell count
%
%   REQUIRES:
%       TopoToolbox 2.x (for GRIDobj resample method)
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

if nargin < 2 || isempty(targetRes)
    targetRes = 500;
end
if nargin < 3
    verbose = true;
end

info = struct();
info.nativeRes = DEM.cellsize;
info.scaleFactor = max(1, round(targetRes / info.nativeRes));
info.coarseRes = info.nativeRes * info.scaleFactor;
info.speedup = info.scaleFactor^2;
info.nativeCells = numel(DEM.Z);

if verbose
    fprintf('  Native resolution:  %.1f m (%d x %d = %.2f M cells)\n', ...
        info.nativeRes, size(DEM.Z, 1), size(DEM.Z, 2), info.nativeCells/1e6);
end

if info.scaleFactor > 1
    DEM_coarse = resample(DEM, info.coarseRes, 'bilinear');
    info.coarseCells = numel(DEM_coarse.Z);
    
    if verbose
        fprintf('  Coarse resolution:  %.1f m (%d x %d = %.2f M cells)\n', ...
            info.coarseRes, size(DEM_coarse.Z, 1), size(DEM_coarse.Z, 2), ...
            info.coarseCells/1e6);
        fprintf('  Speedup factor:     ~%.0fx fewer cells\n', info.speedup);
    end
else
    DEM_coarse = DEM;
    info.coarseCells = info.nativeCells;
    
    if verbose
        fprintf('  Already at or below target resolution. No resampling.\n');
    end
end

end
