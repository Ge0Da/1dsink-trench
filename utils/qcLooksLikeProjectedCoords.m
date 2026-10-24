function tf = qcLooksLikeProjectedCoords(x, y)
%QCLOOKSLIKEPROJECTEDCOORDS Heuristic: coordinates in projected metres.
%
%   Uses the same axis-scale threshold as run_master_figure.m
%   (max |coordinate| > 10000).

x = double(x(:));
y = double(y(:));
if isempty(x) && isempty(y)
    tf = false;
    return;
end
vals = [x; y];
if isempty(vals)
    tf = false;
    return;
end
tf = max(abs(vals), [], 'omitnan') > 10000;
end
