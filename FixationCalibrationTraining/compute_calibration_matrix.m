function [bx, by] = compute_calibration_matrix(tr)
  H_voltage = nan * ones(size(tr));
  H_pixels  = nan * ones(size(tr));
  V_voltage = nan * ones(size(tr));
  V_pixels  = nan * ones(size(tr));
  
  for i = 1:size(tr,2)
    H_voltage(i) = tr(i).H_voltage_median;
    H_pixels(i)  = tr(i).rect_center_X;
    V_voltage(i) = tr(i).V_voltage_median;
    V_pixels(i)  = tr(i).rect_center_Y;
  end

  % add offset 
  H_voltage = [ones(length(H_voltage),1),H_voltage'];
  V_voltage = [ones(length(V_voltage),1),V_voltage'];
  
  % solve
  bx = H_voltage\H_pixels';
  by = V_voltage\V_pixels';  
  
end