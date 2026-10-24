function qcAssert1DSinkCrs(sampleX, sampleY, axisX, axisY)
%QCASSERT1DSINKCRS Fail-closed gate: sample and axis CRS must agree.
%
%   qcAssert1DSinkCrs(sampleX, sampleY, axisX, axisY)
%
%   Errors when geographic sample coordinates would be projected onto a
%   metric trench axis (or the inverse) without prior reprojection.
%   Heuristics match run_master_figure.m.

sampleGeo = qcLooksLikeGeographicCoords(sampleX, sampleY);
sampleProj = qcLooksLikeProjectedCoords(sampleX, sampleY);
axisGeo = qcLooksLikeGeographicCoords(axisX, axisY);
axisProj = qcLooksLikeProjectedCoords(axisX, axisY);

if sampleGeo && axisProj
    error('qc1DSink:crsMismatch', ...
        ['Sample coordinates look geographic (lon/lat) but the trench axis ' ...
         'looks projected (metres). Reproject samples to the DEM CRS before ' ...
         'the 1-D sink path (see reprojectSamplesToGrid).']);
end

if sampleProj && axisGeo
    error('qc1DSink:crsMismatch', ...
        ['Sample coordinates look projected (metres) but the trench axis ' ...
         'looks geographic (lon/lat). Use a projected DEM/axis or supply ' ...
         'sample coordinates in the same CRS.']);
end
end
