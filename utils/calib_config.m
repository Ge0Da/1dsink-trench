function cols = calib_config()
%CALIB_CONFIG  Sample-table columns forming the TI calibration target.
%   The composite Trapping Index is calibrated against the mean normalized
%   content of these heavy metals. Edit this single list to change the target
%   for the whole pipeline (e.g. remove 'HM_Cd_ppm' for the Hg-As-Pb
%   sensitivity test); calibrate_TI_weights.m, geochem_crosscheck.m, and
%   paper_figures.m all read it, so the choice lives in exactly one place.
  cols = {'HM_Hg_ppb', 'HM_As_ppm', 'HM_Cd_ppm', 'HM_Pb_ppm'};
end
