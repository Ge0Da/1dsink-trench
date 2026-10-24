%% [check] qc1DSink — missing DEM gate (fail-closed)
%
%  Run from project root:
%    addpath('utils'); addpath('tests'); check_qc1DSink_missingDem

scriptDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(scriptDir);
addpath(fullfile(projRoot, 'utils'));
addpath(scriptDir);

assertErrorId(@() qcAssert1DSinkDem(''), 'qc1DSink:missingDem');

missingPath = fullfile(projRoot, 'data', '__no_such_dem__.tif');
assertErrorId(@() qcAssert1DSinkDem(missingPath), 'qc1DSink:missingDem');

fixture = fullfile(projRoot, 'data', 'kamchatka_samples.xlsx');
if ~isfile(fixture)
    error('check_qc1DSink_missingDem:fixture', 'Expected fixture: %s', fixture);
end
qcAssert1DSinkDem(fixture);

fprintf('[check] check_qc1DSink_missingDem — PASS\n');
