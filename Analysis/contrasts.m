% Select data files
[filenames, pathname] = uigetfile('*.mat', 'Select data files', 'MultiSelect', 'on');

% Check if the user canceled file selection
if isequal(filenames, 0)
    disp('No files selected. Exiting...');
    return;
end

% Convert filenames to cell array if it's not already
if ~iscell(filenames)
    filenames = {filenames};
end

% Initialize arrays to store contrast and probability data
all_contrast = cell(1, length(filenames));
all_probabilities = cell(1, length(filenames));

% Process each selected data file
for file_idx = 1:length(filenames)
    file_path = fullfile(pathname, filenames{file_idx});
    
    % Load data from file
    loaded_data = load(file_path);
    
    % Check if the loaded data is empty or doesn't contain expected fields
    if isempty(loaded_data) || ~isfield(loaded_data, 'expt_info') || ~isfield(loaded_data.expt_info, 'trial_records')
        disp(['Error loading data from file: ', file_path]);
        continue;  % Skip to the next file
    end
    
    trial_records = loaded_data.expt_info.trial_records;
    
    % Check if trial_records is empty
    if isempty(trial_records)
        disp(['No trial records found in file: ', file_path]);
        continue;  % Skip to the next file
    end
    
    % Extract contrast values and probabilities
    contrast = ones(1, length(trial_records)) * nan;
    for i = 1:length(trial_records)
        contrast(i) = trial_records(i).contrast;
    end
    contrast = unique(contrast);
    
    p = ones(1, length(contrast)) * nan;
    for c = 1:length(contrast)
        hitmiss = [];
        for i = 1:length(trial_records)
            if trial_records(i).contrast == contrast(c)
                if strcmp(trial_records(i).trial_error, 'hit')
                    hitmiss = [hitmiss, 1];
                elseif strcmp(trial_records(i).trial_error, 'wrong_target')
                    hitmiss = [hitmiss, 0];
                end
            end
        end
        p(c) = sum(hitmiss) / length(hitmiss);
    end
    
    % Store contrast and probability data for this file
    all_contrast{file_idx} = contrast;
    all_probabilities{file_idx} = p;
end

% Plot contrast vs. probability for each data file
figure;
hold on;
for file_idx = 1:length(filenames)
    plot(all_contrast{file_idx}, all_probabilities{file_idx}, 'o-');
end
hold off;

xlabel('Contrast');
ylabel('Probability');
title('Contrast vs. Probability');

% Add legend
legend(filenames, 'Interpreter', 'none');

