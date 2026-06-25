function [xProj, yProj, crsInfo] = reprojectSamplesToGrid(lon, lat, demFile)
%REPROJECTSAMPLESTOGRID  Reproject lon/lat coordinates into a DEM's CRS.
%
%   [xProj, yProj] = reprojectSamplesToGrid(lon, lat, demFile)
%
%   Reads the coordinate reference system directly from a GeoTIFF file
%   and forward-projects geographic coordinates (longitude, latitude) into
%   the same projected coordinate system the DEM uses.
%
%   This function exists because:
%     - Sample files typically store positions as longitude/latitude (degrees)
%     - DEMs loaded by GRIDobj/TopoToolbox may use any projected CRS
%       (UTM, Mercator, custom Transverse Mercator, etc.)
%     - Mixing degrees and metres produces wrong results silently
%
%   Three fallback paths are tried, in order of robustness:
%     1. georasterinfo  → projcrs object  → projfwd   (R2021b+, any CRS)
%     2. geotiffinfo    → EPSG code       → projcrs   (standard EPSG only)
%     3. geotiffinfo    → geotiff2mstruct → projfwd   (any CRS)
%
%   INPUTS:
%     lon      - Longitude values (degrees East), Nx1 or 1xN
%     lat      - Latitude values  (degrees North), Nx1 or 1xN
%     demFile  - Path to the GeoTIFF file whose CRS to match
%
%   OUTPUTS:
%     xProj    - Projected x-coordinates (metres), same shape as lon
%     yProj    - Projected y-coordinates (metres), same shape as lat
%     crsInfo  - String describing which method succeeded
%
%   ERRORS:
%     Throws an error if all three paths fail (missing Mapping Toolbox,
%     corrupt GeoTIFF, or unreadable CRS metadata).
%
%   EXAMPLE:
%     % In a driver script, after loading samples and DEM:
%     [samples.X, samples.Y] = reprojectSamplesToGrid( ...
%         samples.Longitude, samples.Latitude, 'kamchatka_gmrt.tif');
%
%   REQUIRES: Mapping Toolbox
%
%   See also: projfwd, georasterinfo, geotiffinfo, geotiff2mstruct, projcrs
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2026
%   Version: 1.0

% -------------------------------------------------------------------------
%  Input validation
% -------------------------------------------------------------------------
arguments
    lon  (:,1) double
    lat  (:,1) double
    demFile (1,:) char
end

if ~isfile(demFile)
    error('reprojectSamplesToGrid:fileNotFound', ...
        'DEM file not found: %s', demFile);
end

% -------------------------------------------------------------------------
%  Path 1: georasterinfo → projcrs → projfwd  (most robust, R2021b+)
%
%  georasterinfo reads the FULL projection definition from the GeoTIFF
%  (WKT string), not just an EPSG code.  This handles user-defined
%  projections (PCS = 32767) that have no standard EPSG number —
%  common in GMRT and GMT-produced grids.
% -------------------------------------------------------------------------
try
    info   = georasterinfo(demFile);
    demCRS = info.CoordinateReferenceSystem;
    [xProj, yProj] = projfwd(demCRS, lat, lon);
    crsInfo = 'georasterinfo → projfwd';
    return
catch
end

% -------------------------------------------------------------------------
%  Path 2: geotiffinfo → EPSG code → projcrs → projfwd
%
%  Works when the GeoTIFF stores a standard EPSG projected-CRS code
%  (e.g. 32655 for UTM zone 55N).  Fails silently for user-defined CRS.
%  This path mirrors what ModuleA_Main_Driver does via demInfo.epsgCode.
% -------------------------------------------------------------------------
try
    tiffInfo = geotiffinfo(demFile);
    if isfield(tiffInfo, 'GeoTIFFCodes') && isfield(tiffInfo.GeoTIFFCodes, 'PCS')
        pcs = tiffInfo.GeoTIFFCodes.PCS;
        if pcs > 0 && pcs ~= 32767
            crs = projcrs(pcs);
            [xProj, yProj] = projfwd(crs, lat, lon);
            crsInfo = sprintf('EPSG:%d → projfwd', pcs);
            return
        end
    end
catch
end

% -------------------------------------------------------------------------
%  Path 3: geotiffinfo → geotiff2mstruct → projfwd
%
%  Builds a MATLAB map projection structure from the full GeoTIFF
%  metadata.  Works even for user-defined projections on older MATLAB
%  versions that lack georasterinfo.
% -------------------------------------------------------------------------
try
    if ~exist('tiffInfo', 'var')
        tiffInfo = geotiffinfo(demFile);
    end
    mstruct = geotiff2mstruct(tiffInfo);
    [xProj, yProj] = projfwd(mstruct, lat, lon);
    crsInfo = 'geotiff2mstruct → projfwd';
    return
catch
end

% -------------------------------------------------------------------------
%  All paths failed
% -------------------------------------------------------------------------
error('reprojectSamplesToGrid:allPathsFailed', ...
    ['Could not read CRS from ''%s''.\n' ...
     'Ensure the Mapping Toolbox is installed and the file\n' ...
     'is a valid GeoTIFF with embedded projection metadata.\n' ...
     'Alternatively, supply pre-projected sample coordinates.'], ...
    demFile);

end
