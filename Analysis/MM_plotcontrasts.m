% Select data files
[filename, pathname] = uigetfile('*.mat', 'Select data files', 'MultiSelect', 'on');

% Check if the user canceled file selection
if isequal(filename, 0)
    disp('No files selected. Exiting...');
    return;
end

% Convert filename to cell array if it's not already
if ~iscell(filename)
    filename = {filename};
end

% Concatenate the path and filenames
data_files = fullfile(pathname, filename);

% Initialize arrays to store contrast and probability data
all_contrast = [];
all_probabilities = [];

% Process each selected data file
for file_idx = 1:length(data_files)
    file_path = data_files{file_idx};
   
    % Load data from file
    loaded_data = load(file_path); 

    trial_records = loaded_data.expt_info.trial_records;

    % Extract contrast values from trial records
    contrast_values = [];
    for i = 1:length(trial_records)
        contrast_values(i) = trial_records(i).contrast;
    end

    % Extract hit/miss information and calculate probability for each contrast
    unique_contrasts = unique(contrast_values);
    probabilities = zeros(size(unique_contrasts));
    for c = 1:length(unique_contrasts)
        hit_miss = [];
        for i = 1:length(trial_records)
            if trial_records(i).contrast == unique_contrasts(c)
                if strcmp(trial_records(i).trial_error, 'hit')
                    hit_miss = [hit_miss, 1];
                elseif strcmp(trial_records(i).trial_error, 'wrong target')
                    hit_miss = [hit_miss, 0];
                end
            end
        end
        probabilities(c) = sum(hit_miss) / length(hit_miss);
    end

    % Store contrast and probability data for this file
    all_contrast = [all_contrast, unique_contrasts];
    all_probabilities = [all_probabilities, probabilities];
end

% Plot contrast vs probability
figure;
plot(all_contrast, all_probabilities, 'o');
xlabel('Contrast');
ylabel('Probability');

% Use Weibull function to fit the data and extract contrast at 75% performance
[params, ~] = fminsearch(@(params) weibull_fit(params, all_contrast, all_probabilities), [0.5, 1]);
contrast_75 = weibull_inverse(params, 0.75);

% Display the contrast at 75% performance
disp(['Contrast at 75% performance: ', num2str(contrast_75)]);

% Weibull functions for fitting
function error = weibull_fit(params, contrast, probabilities)
    contrast_75 = weibull_inverse(params, 0.75);
    error = sum((probabilities - weibull_function(params, contrast)).^2);
end

function y = weibull_function(params, x)
    y = 1 - exp(-(x / params(1)).^params(2));
end

function x = weibull_inverse(params, y)
    x = params(1) * (-log(1 - y)).^(1 / params(2));
end
