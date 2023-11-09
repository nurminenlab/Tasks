function [PC,response,running_PC] = percentage_correct(trial_records)
  
  # for moving average
  n_trials = 10;
  response = nan*ones(1,length(trial_records));
  
  for i = 1:length(trial_records)
    if strcmp(trial_records(i).trial_error,'hit')
      response(i) = 1;
    elseif strcmp(trial_records(i).trial_error,'miss')
      response(i) = 0;
    else
      # do nothing
    end  
  endfor
  response = response(~isnan(response));
  PC = sum(response) / length(response);
  B = (1/n_trials) * ones(n_trials,1);
  running_PC = filter(B,1,response);
