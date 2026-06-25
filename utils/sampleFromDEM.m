function z = sampleFromDEM(x, y, DEM)
%SAMPLEFROMDEM Sample elevation values from DEM at given coordinates
%
%   z = sampleFromDEM(x, y, DEM)
%
%   Part of S2S ModuleA - Source-to-Sink Geomorphometric Analysis
%
%   Simple wrapper around interp() for clarity and consistency.
%
%   INPUTS:
%       x, y - Coordinate vectors
%       DEM  - GRIDobj to sample from
%
%   OUTPUTS:
%       z    - Sampled elevation values
%
%   REQUIRES:
%       TopoToolbox 2.x
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

z = interp(DEM, x(:), y(:));

end
