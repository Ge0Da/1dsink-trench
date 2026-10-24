function qcRun1DSinkPathGates(stage, ctx)
%QCRUN1DSINKPATHGATES Run fail-closed QC gates for the 1-D sink path.
%
%   qcRun1DSinkPathGates('dem', struct('demPath', demPath))
%   qcRun1DSinkPathGates('axis', struct('trenchAxis', trenchAxis))
%   qcRun1DSinkPathGates('samples', struct( ...
%       'sampleX', sp_x, 'sampleY', sp_y, ...
%       'axisX', trenchAxis.x, 'axisY', trenchAxis.y))
%
%   STAGES (fail-closed — any failure throws):
%     'dem'     — DEM file present
%     'axis'    — no excessive NaN runs on trench axis x/y/z
%     'samples' — sample vs axis CRS consistency

arguments
    stage (1,:) char {mustBeMember(stage, {'dem', 'axis', 'samples'})}
    ctx struct
end

switch stage
    case 'dem'
        if ~isfield(ctx, 'demPath')
            error('qc1DSink:badContext', 'dem stage requires ctx.demPath.');
        end
        qcAssert1DSinkDem(ctx.demPath);

    case 'axis'
        if ~isfield(ctx, 'trenchAxis')
            error('qc1DSink:badContext', 'axis stage requires ctx.trenchAxis.');
        end
        qcAssert1DSinkNanRuns(ctx.trenchAxis);

    case 'samples'
        need = {'sampleX', 'sampleY', 'axisX', 'axisY'};
        for k = 1:numel(need)
            if ~isfield(ctx, need{k})
                error('qc1DSink:badContext', ...
                    'samples stage requires ctx.%s.', need{k});
            end
        end
        qcAssert1DSinkCrs(ctx.sampleX, ctx.sampleY, ctx.axisX, ctx.axisY);
end
end
