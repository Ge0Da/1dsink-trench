%% [check] qc1DSink — CRS mismatch gate (fail-closed)
%
%  Run from project root:
%    addpath('utils'); addpath('tests'); check_qc1DSink_crsMismatch

scriptDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(scriptDir);
addpath(fullfile(projRoot, 'utils'));
addpath(scriptDir);

lon = [155; 160];
lat = [45; 50];
axisX = 5e5 + (0:4)' * 1000;
axisY = 5.5e6 + (0:4)' * 1000;

assertErrorId(@() qcAssert1DSinkCrs(lon, lat, axisX, axisY), 'qc1DSink:crsMismatch');

projX = 5e5 + (0:2)' * 100;
projY = 5.5e6 + (0:2)' * 100;
qcAssert1DSinkCrs(projX, projY, axisX(1:3), axisY(1:3));

geoAxisX = lon;
geoAxisY = lat;
assertErrorId(@() qcAssert1DSinkCrs(projX, projY, geoAxisX, geoAxisY), ...
    'qc1DSink:crsMismatch');

fprintf('[check] check_qc1DSink_crsMismatch — PASS\n');
