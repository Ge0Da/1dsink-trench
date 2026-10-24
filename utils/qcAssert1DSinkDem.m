function qcAssert1DSinkDem(demPath)
%QCASSERT1DSINKDEM Fail-closed gate: DEM path must exist on disk.
%
%   qcAssert1DSinkDem(demPath)
%
%   Part of 1DSink Along-Trench Analysis Pipeline — ZONE1 methods hygiene.

arguments
    demPath (1,:) {mustBeText}
end

demPath = strtrim(char(demPath));
if isempty(demPath)
    error('qc1DSink:missingDem', 'DEM path is empty.');
end

if ~isfile(demPath)
    error('qc1DSink:missingDem', 'DEM file not found: %s', demPath);
end
end
