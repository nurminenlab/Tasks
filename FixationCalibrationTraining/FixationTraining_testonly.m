% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;
clear;

% user defined parameters
scaler = 0.5;
stimulus_duration = 1;
isi_duration = 1;
ms = 10;

[FNAME, FPATH, FLTIDX] = uigetfile();
load([FPATH,FNAME]);
[bx,by] = compute_calibration_matrix(tr);

% set-up Datapixx
Datapixx('Open');
adcRate = 1e3;
nAdcSamples = floor((stimulus_duration + isi_duration) * adcRate);
minStreamFrames = floor(adcRate / 50);

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if available
eyeTrack_screenNumber = min(screens);
stimulus_screenNumber = max(screens);

% Define black and white
white = WhiteIndex(stimulus_screenNumber);
black = BlackIndex(stimulus_screenNumber);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[stimulus_window, stimulus_windowRect] = PsychImaging('OpenWindow', stimulus_screenNumber, grey);
[eyeTrack_window, eyeTrack_windowRect] = PsychImaging('OpenWindow', eyeTrack_screenNumber, grey);

try 
  [x,y,buttons,focus] = GetMouse(eyeTrack_window);
catch
  sca();
end

sca();
