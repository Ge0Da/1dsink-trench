%% run_check_tests — leaf-node [check] QC tests for the 1-D sink path
%
%  addpath('utils'); addpath('tests'); run_check_tests

scriptDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(scriptDir);
addpath(fullfile(projRoot, 'utils'));
addpath(scriptDir);

checks = {
    'check_qc1DSink_missingDem'
    'check_qc1DSink_nanRuns'
    'check_qc1DSink_crsMismatch'
};

fprintf('\n[check] Running %d QC gate tests...\n\n', numel(checks));
for k = 1:numel(checks)
    run(checks{k});
end
fprintf('\n[check] All QC gate tests passed.\n\n');
