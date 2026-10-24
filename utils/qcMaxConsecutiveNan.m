function nMax = qcMaxConsecutiveNan(v)
%QCMAXCONSECUTIVENAN Longest run of consecutive NaNs in a vector.

v = v(:);
if isempty(v)
    nMax = 0;
    return;
end

isNan = isnan(v);
if ~any(isNan)
    nMax = 0;
    return;
end

d = diff([0; isNan; 0]);
runStarts = find(d == 1);
runEnds = find(d == -1) - 1;
nMax = max(runEnds - runStarts + 1);
end
