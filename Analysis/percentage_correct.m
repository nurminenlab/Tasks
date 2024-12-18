function [PC,response,running_PC] = percentage_correct(expt_info)

  trial_records = expt_info.trial_records;  
  # for moving average
  n_trials = 20;
  response = nan*ones(1,length(trial_records));
  
  for i = 1:length(trial_records)
    if strcmp(trial_records(i).trial_error,'hit') || strcmp(trial_records(i).trial_error,'wrong_target')
      response(i) = 1;
    else
      response(i) = 0;    
      # do nothing
    end  
  endfor
  response = response(~isnan(response));
  PC = sum(response) / length(response);
  B = (1/n_trials) * ones(n_trials,1);
  running_PC = filter(B,1,response);
