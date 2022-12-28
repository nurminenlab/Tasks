% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;
clear;

PsychDefaultSetup(2);
PsychDebugWindowConfiguration;
screens = Screen('Screens');
stimulus_screenNumber = max(screens);

white = WhiteIndex(stimulus_screenNumber);
black = BlackIndex(stimulus_screenNumber);
grey = 128;
inc = white - grey;

stimulus_window = Screen('OpenWindow',stimulus_screenNumber,128);
[screenXpixels, screenYpixels] = Screen('WindowSize',stimulus_window);

gridSize = 256;
parameters.lowCut_S  = 0.01;
parameters.highCut_S = 0.05;
parameters.orientation_low_S  = 0;
parameters.orientation_high_S = 90;
parameters.plateauPixels_S = 200;
parameters.edgePixels_S = 20;
parameters.contrast_S = 0.5;

parameters.lowCut_C  = 0.2;
parameters.highCut_C = 0.25;
parameters.orientation_low_C  = 135-25;
parameters.orientation_high_C = 135+25;
parameters.orientations_C = [0,90];
parameters.plateauPixels_C = 100;
parameters.edgePixels_C = 1;
parameters.contrast_C = 0.5;


#[grText] = generate_filteredNoise_CS(gridSize,stimulus_window,parameters);
[grText] = generate_filteredNoise(gridSize,stimulus_window,parameters);


rect = [0 0 gridSize gridSize];
rect = CenterRectOnPoint(rect, screenXpixels/2, screenYpixels/2);
rects  = nan * ones(4,4);
thetas = deg2rad([-45,45,135,225]);
R = 200;
[Xoff, Yoff] = pol2cart(thetas,R);
for i = 1:length(Xoff)
  rects(:,i) = round(CenterRectOnPoint(rect, screenXpixels/2+Xoff(i), screenYpixels/2+Yoff(i)));
end


Screen('DrawTextures', stimulus_window, grText, [], rects);
Screen('Flip', stimulus_window);
KbWait();
sca;