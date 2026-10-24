function qcAssert1DSinkNanRuns(axis, opts)
%QCASSERT1DSINKNANRUNS Fail-closed gate: limit consecutive NaNs on axis.
%
%   qcAssert1DSinkNanRuns(axis)
%   qcAssert1DSinkNanRuns(axis, 'maxConsecutiveNan', 1)
%
%   axis must contain fields x, y, z (from extractTrenchAxis).

arguments
    axis struct
    opts.maxConsecutiveNan (1,1) double {mustBeNonnegative} = 1
end

requiredFields = {'x', 'y', 'z'};
for k = 1:numel(requiredFields)
    f = requiredFields{k};
    if ~isfield(axis, f)
        error('qc1DSink:badAxis', 'axis missing field ''%s''.', f);
    end
end

fieldLabels = {'x', 'y', 'z'};
for k = 1:numel(fieldLabels)
    f = fieldLabels{k};
    nRun = qcMaxConsecutiveNan(axis.(f));
    if nRun > opts.maxConsecutiveNan
        error('qc1DSink:nanRun', ...
            ['Along-trench axis field ''%s'' has %d consecutive NaN values ' ...
             '(limit %d). Aborting 1-D sink path.'], ...
            f, nRun, opts.maxConsecutiveNan);
    end
end
end
