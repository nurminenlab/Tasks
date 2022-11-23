function [bx, by, HV, VV, HP,VP] = compute_calibration_matrix(tr, plot_yes)
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
  
  HV = H_voltage(:,2);
  VV = V_voltage(:,2);  
  HP = H_pixels;
  VP = V_pixels; 
 
  
  if plot_yes == 1
    hf1 = figure(1);
    hold on     
    XX = linspace(min(HV), max(HV),4);
    plot(HV, HP,'ko')
    plot(XX,bx(1) + bx(2)*XX,'k-')
    
    hf2 = figure(2);
    hold on 
    YY = linspace(min(VV), max(VV),4);
    plot(VV, VP,'ko')        
    plot(YY,by(1) + by(2)*YY,'k-')
    
    KbStrokeWait();
    close(hf1);
    close(hf2);    
  endif
end