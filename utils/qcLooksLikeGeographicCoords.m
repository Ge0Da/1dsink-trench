function tf = qcLooksLikeGeographicCoords(x, y)
%QCLOOKSLIKEGEOGRAPHICCOORDS Heuristic: coordinates in lon/lat degrees.
%
%   Matches the sample-vs-axis CRS check used in run_master_figure.m
%   (all |x| <= 360 and all |y| <= 90).

x = double(x(:));
y = double(y(:));
if isempty(x) || isempty(y) || numel(x) ~= numel(y)
    tf = false;
    return;
end
tf = all(abs(x) <= 360) && all(abs(y) <= 90);
end
