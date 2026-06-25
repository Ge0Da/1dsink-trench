function hm = targetMean(T, colNames)
%TARGETMEAN  Min-max normalized mean of the TI calibration-target metals.
%   hm = TARGETMEAN(T, colNames) returns, for each row of sample table T, the
%   mean over the metals listed in calib_config() after scaling each metal to
%   [0, 1]. Rows missing every target metal return NaN. This is the single
%   construction of the calibration target shared across the pipeline.
  cols = calib_config();
  M = nan(height(T), numel(cols));
  for i = 1:numel(cols)
    v = readCol(T, colNames, cols{i});
    lo = min(v, [], 'omitnan');
    hi = max(v, [], 'omitnan');
    M(:, i) = (v - lo) ./ (hi - lo + eps);
  end
  hm = mean(M, 2, 'omitnan');
end

function v = readCol(Tbl, names, colname)
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
