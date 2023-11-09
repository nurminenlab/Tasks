addpath('/home/vpixx/Tasks/Functions/');
addpath('/home/vpixx/Tasks/Analysis/');

dataroot = '/home/vpixx/TrialRecords/Tully_trialrecords/';
# first reasonable Texture task performance
load([dataroot,'Tully-TextureTask-trial_records-16-Mar-2023.mat']);

gridSize = 256;
parameters.lowCut_S  = 0.01;
parameters.highCut_S = 0.05;
parameters.orientation_low_S  = 135-45;
parameters.orientation_high_S = 135+45;
parameters.plateauPixels_S = 200;
parameters.edgePixels_S = 20;
parameters.contrast_S = 0.5;

parameters.lowCut_C  = 0.2;
parameters.highCut_C = 0.25;
parameters.orientation_low_C  = 45-25;
parameters.orientation_high_C = 45+25;
#parameters.orientations_C = [0,90];
parameters.plateauPixels_C = 100;
parameters.edgePixels_C = 1;
parameters.contrast_C = 1;


[grayScaleImageMatrix_CS, grayScaleImageMatrix_S] = generate_filteredNoise_CS_illustrations(gridSize,parameters);

LUT = [1 1 1]'*linspace(0,1,256);

figure();
image(grayScaleImageMatrix_CS); axis square off; colormap(LUT');
figure();
image(grayScaleImageMatrix_S); axis square off; colormap(LUT');

[PC,response,running_PC] = percentage_correct(trial_records);

figure();
hold on 
plot(1:length(response), running_PC,'r-','LineWidth',2,'Color',[1 0 0])
plot([1,length(response)], [0.25,0.25],'k--','LineWidth',2,'Color',[0 0 0])
plot([1,length(response)], [PC,PC],'b-','LineWidth',4,'Color',[0 0 1])
box off;
set(gca,'tickdir','out')
xlabel('Trial number')
ylabel('Proportion correct')

# load fixation data for large faces
YBC = dir([dataroot,'Tully_trial_records*']);
on_target_time = [];
for i = 1:length(YBC)
  load([dataroot,YBC(i).name]);
  for j = 1:length(trial_records)
    on_target_time = [on_target_time,trial_records(j).on_target_time];      
  end  
end

figure();
on_target_time(on_target_time < 0.1) = [];
hist(on_target_time,0.1:0.1:5)
box off
set(gca,'tickdir','out');
xlabel('Fixation duration (s)')
ylabel('N trials');








