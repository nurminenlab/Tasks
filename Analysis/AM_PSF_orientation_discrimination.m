trial_records=expt_info.trial_records;
contrast=ones(1,length(trial_records))*nan;
for i = 1:length(trial_records)
  contrast(i) = trial_records(i).contrast;  
end
contrast = unique(contrast);
p = ones(1,length(contrast))*nan;
for c = 1:length(contrast)
  hitmiss = [];
  for i = 1:length(trial_records)
    if trial_records(i).contrast == contrast(c)
      if strcmp(trial_records(i).trial_error,'hit')
        hitmiss = [hitmiss,1];
      elseif strcmp(trial_records(i).trial_error,'wrong_target')
        hitmiss = [hitmiss,0];
      else
        # do nothing
      endif
    endif
  endfor  
  p(c) = sum(hitmiss)/length(hitmiss);
end
plot(contrast,p,'o-');




 