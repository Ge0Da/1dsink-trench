%% [check] qc1DSink — NaN run gate (fail-closed)
%
%  Run from project root:
%    addpath('utils'); addpath('tests'); check_qc1DSink_nanRuns

scriptDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(scriptDir);
addpath(fullfile(projRoot, 'utils'));
addpath(scriptDir);

axisOk = struct('x', (1:5)', 'y', (1:5)', 'z', -(1:5)');
qcAssert1DSinkNanRuns(axisOk);

axisBad = struct('x', [1; 2; nan; nan; 5], 'y', (1:5)', 'z', -(1:5)');
assertErrorId(@() qcAssert1DSinkNanRuns(axisBad), 'qc1DSink:nanRun');

assert(qcMaxConsecutiveNan([1; nan; 3]) == 1);
assert(qcMaxConsecutiveNan([1; nan; nan; 4]) == 2);

fprintf('[check] check_qc1DSink_nanRuns — PASS\n');
