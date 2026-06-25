% geochem_crosscheck.m
% =========================================================================
%  Cross-correlate the Trapping Index (TI) with sediment geochemistry
% =========================================================================
%
%  Validates the morphometric Trapping Index against independent contaminant
%  concentrations in Kuril-Kamchatka Trench surface sediments. Reports:
%    1. Spearman correlations of TI vs each contaminant class and per element
%       (headline target: four-element heavy-metal mean Hg, As, Cd, Pb)
%    2. Multiple-comparison correction (Benjamini-Hochberg) on the secondary
%       per-element and per-class tests
%    3. Partial correlations controlling for along-axis position (source vs
%       morphology) and for sediment mud fraction and TOC (grain-size proxy)
%
%  Inputs:
%    output/1DSink_MasterFigure_results.mat   TI and projected sample positions
%    data/kamchatka_samples.xlsx              contaminant concentrations
%    data/sediment_texture_TOC.csv            per-station mud% and TOC
%    output/1DSink_projected_samples.csv      along-axis and lateral positions
%
%  Output: console table + output/TI_contaminant_correlations.csv
%
%  Author: David Fernandez-Blanco
% =========================================================================

if ~exist('tTotal', 'var'), clear; clc; end
projRoot = fileparts(fileparts(mfilename('fullpath')));
if isempty(projRoot), projRoot = fileparts(pwd); end
addpath(fullfile(projRoot, 'utils'));
cd(projRoot);

%% 1. Load TI and interpolate to sample positions -------------------------
ld = load('output/1DSink_MasterFigure_results.mat');
mr = ld.masterResults;
TI_s_km = ld.bands.s_km(:);
TI_val  = mr.metrics.TI_full(:);
vTI     = ~isnan(TI_val) & ~isnan(TI_s_km);

ps          = ld.projectedSamples;
sample_s_km = ps.s_km(:);
nSamples    = numel(sample_s_km);
TI = interp1(TI_s_km(vTI), TI_val(vTI), sample_s_km, 'linear', NaN);

% Station names (for joining the auxiliary tables).
if isfield(ps, 'name'), rawNames = ps.name; else, rawNames = ps.Name; end
stn = cell(nSamples, 1);
for ii = 1:nSamples
    v = rawNames(ii); if iscell(v), v = v{1}; end
    stn{ii} = strtrim(char(v));
end

%% 2. Load contaminants ----------------------------------------------------
T = readtable('data/kamchatka_samples.xlsx', 'VariableNamingRule', 'preserve');
cn = T.Properties.VariableNames;

Total_OCP = sumOrNaN(table2array(T(:, startsWith(cn, 'OCP_'))));
Total_PAH = sumOrNaN(table2array(T(:, startsWith(cn, 'PAH_'))));
hydAll    = startsWith(cn, 'HYD_');
excl      = contains(cn, 'D50_IS') | contains(cn, '3_methyl') | ...
            contains(cn, 'Pristan') | contains(cn, 'Phytane') | ...
            contains(cn, 'hexadecene') | contains(cn, 'octadecene');
Total_HYD = sumOrNaN(table2array(T(:, hydAll & ~excl)));   % strict n-alkanes C14-C29

Hg = getCol(T, cn, 'HM_Hg_ppb');
As = getCol(T, cn, 'HM_As_ppm');
Cd = getCol(T, cn, 'HM_Cd_ppm');
Pb = getCol(T, cn, 'HM_Pb_ppm');

HM = targetMean(T, cn);   % calibration target (calib_config; default Hg, As, Cd, Pb)

%% 3. Spearman correlations + multiple-comparison correction --------------
labels  = {'HM mean (Hg,As,Cd,Pb)', 'Hg', 'As', 'Cd', 'Pb', ...
           'PAH', 'OCP', 'n-alkanes'};
vecs    = {HM, Hg, As, Cd, Pb, Total_PAH, Total_OCP, Total_HYD};
rho = nan(numel(vecs), 1); pval = rho; nObs = rho;
for i = 1:numel(vecs)
    [rho(i), pval(i), nObs(i)] = spear(TI, vecs{i});
end

fprintf('\nTI vs contaminant (Spearman)\n%s\n', repmat('-', 1, 52));
fprintf('%-20s %8s %8s %5s\n', 'Quantity', 'rho', 'p', 'n');
for i = 1:numel(vecs)
    fprintf('%-20s %+8.3f %8.4f %5d\n', labels{i}, rho(i), pval(i), nObs(i));
end

% Secondary family (per-element + organics; the four-element mean is the primary
% hypothesis and is not corrected). Benjamini-Hochberg false-discovery rate.
secIdx = 2:numel(vecs);
qBH = bhFDR(pval(secIdx));
fprintf('\nMultiple-comparison correction (BH-FDR), secondary tests\n%s\n', repmat('-', 1, 52));
fprintf('%-20s %8s %8s %10s\n', 'Test', 'rho', 'p', 'q (BH)');
for k = 1:numel(secIdx)
    i = secIdx(k);
    fprintf('%-20s %+8.3f %8.4f %10.4f\n', labels{i}, rho(i), pval(i), qBH(k));
end

%% 4. Partial correlations -- along-axis position (source vs morphology) ---
PS = readtable('output/1DSink_projected_samples.csv', 'VariableNamingRule', 'preserve');
[along, lateral] = deal(nan(nSamples, 1));
for ii = 1:nSamples
    h = find(strcmp(strtrim(string(PS.Name)), stn{ii}), 1);
    if ~isempty(h), along(ii) = PS.AlongAxis_km(h); lateral(ii) = PS.DistToAxis_km(h); end
end

fprintf('\nPartial Spearman -- control along-axis position\n%s\n', repmat('-', 1, 64));
fprintf('%-16s %18s %18s %18s\n', 'Pair', 'zero-order', '| position', '| position+lateral');
for pr = {'HM mean', HM; 'Hg', Hg; 'As', As; 'Pb', Pb}'
    y = pr{2};
    v0 = ~isnan(TI) & ~isnan(y);
    [r0, p0] = corr(TI(v0), y(v0), 'type', 'Spearman');
    v1 = v0 & ~isnan(along);
    [r1, p1] = partialcorr(TI(v1), y(v1), along(v1), 'type', 'Spearman');
    v2 = v1 & ~isnan(lateral);
    [r2, p2] = partialcorr(TI(v2), y(v2), [along(v2) lateral(v2)], 'type', 'Spearman');
    fprintf('%-16s %+8.3f (%.3f) %+8.3f (%.3f) %+8.3f (%.3f)\n', pr{1}, r0, p0, r1, p1, r2, p2);
end

%% 5. Partial correlations -- grain size and TOC (grain-size proxy) --------
if isfile('data/sediment_texture_TOC.csv')
    GS = readtable('data/sediment_texture_TOC.csv', 'VariableNamingRule', 'preserve');
    [mud, toc] = deal(nan(nSamples, 1));
    for ii = 1:nSamples
        h = find(strcmp(strtrim(string(GS.Site)), stn{ii}), 1);
        if ~isempty(h), mud(ii) = GS.Mud_pct(h); toc(ii) = GS.TOC_pct(h); end
    end
    fprintf('\nPartial Spearman -- control mud%% and TOC\n%s\n', repmat('-', 1, 70));
    fprintf('%-16s %14s %14s %14s %14s  n\n', 'Pair', 'zero-order', '| mud%', '| TOC', '| mud%+TOC');
    for pr = {'HM mean', HM; 'Hg', Hg; 'As', As; 'Pb', Pb}'
        y = pr{2};
        vb = ~isnan(TI) & ~isnan(y) & ~isnan(mud) & ~isnan(toc);
        [r0, ~]  = corr(TI(vb), y(vb), 'type', 'Spearman');
        [rm, pm] = partialcorr(TI(vb), y(vb), mud(vb), 'type', 'Spearman');
        [rt, pt] = partialcorr(TI(vb), y(vb), toc(vb), 'type', 'Spearman');
        [rbth, ~] = partialcorr(TI(vb), y(vb), [mud(vb) toc(vb)], 'type', 'Spearman');
        fprintf('%-16s %+8.3f %+8.3f (%.3f) %+8.3f (%.3f) %+8.3f  %d\n', ...
            pr{1}, r0, rm, pm, rt, pt, rbth, sum(vb));
    end
end

%% 6. Write a compact correlation table -----------------------------------
out = table(labels(:), rho, pval, nObs, 'VariableNames', {'Quantity', 'rho', 'p', 'n'});
writetable(out, 'output/TI_contaminant_correlations.csv');
fprintf('Saved output/TI_contaminant_correlations.csv\n\n');

%% ------------------------------------------------------------------------
function v = sumOrNaN(M)
  if isempty(M), v = nan(size(M, 1), 1); return; end
  v = sum(M, 2, 'omitnan');
  v(all(isnan(M), 2)) = NaN;
end

function v = getCol(Tbl, names, colname)
  if any(strcmp(names, colname))
    raw = Tbl.(colname);
    if iscell(raw)
      v = nan(numel(raw), 1);
      for jj = 1:numel(raw)
        if isnumeric(raw{jj}) && isscalar(raw{jj}), v(jj) = double(raw{jj});
        elseif ischar(raw{jj}) || isstring(raw{jj}), v(jj) = str2double(raw{jj}); end
      end
    else
      v = double(raw);
    end
  else
    v = nan(height(Tbl), 1);
  end
end

function [rho, p, n] = spear(x, y)
  valid = ~isnan(x) & ~isnan(y); n = sum(valid);
  if n < 4, rho = NaN; p = NaN; return; end
  [rho, p] = corr(x(valid), y(valid), 'type', 'Spearman');
end

function q = bhFDR(p)
  p = p(:); m = numel(p); [ps, ord] = sort(p, 'ascend');
  qs = ps .* m ./ (1:m)';
  for kk = m-1:-1:1, qs(kk) = min(qs(kk), qs(kk+1)); end
  qs = min(qs, 1); q = nan(m, 1); q(ord) = qs;
end
