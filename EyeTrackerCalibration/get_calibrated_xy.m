function [cx,cy] = get_calibrated_xy(tr,bx,by)
  
  cx = nan*ones(length(tr));
  cy = nan*ones(length(tr));
  
  for i = 1:length(tr)
    cx(i) = tr(i).H_voltage_median*bx(2) + bx(1);
    cy(i) = tr(i).V_voltage_median*by(2) + by(1);
  end
  
endfunction
