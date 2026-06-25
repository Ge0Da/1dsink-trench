function [FD, A, DEM_filled] = computeFlowRouting(DEM, verbose)
%COMPUTEFLOWROUTING Standard flow routing: fill → flow direction → accumulation
%
%   [FD, A, DEM_filled] = computeFlowRouting(DEM, verbose)
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   INPUTS:
%       DEM     - GRIDobj (typically at coarse resolution for speed)
%       verbose - Print progress messages [default: true]
%
%   OUTPUTS:
%       FD         - FLOWobj
%       A          - Flow accumulation (GRIDobj)
%       DEM_filled - Sink-filled DEM
%
%   REQUIRES:
%       TopoToolbox 2.x
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 1.0

if nargin < 2
    verbose = true;
end

if verbose
    fprintf('  Filling sinks...\n');
end
DEM_filled = fillsinks(DEM);

if verbose
    fprintf('  Computing flow direction...\n');
end
FD = FLOWobj(DEM_filled, 'preprocess', 'carve');

if verbose
    fprintf('  Computing flow accumulation...\n');
end
A = flowacc(FD);

end
