function [trans_mtrx,scale_mtrx] = set_tracker_GainOffset(debug_on)

# load fixation data
[FNAME, FPATH, FLTIDX] = uigetfile([],'Select gaze data');
load([FPATH,FNAME]);
# load rough calibration
[FNAME, FPATH, FLTIDX] = uigetfile('Select manual calibration');
load([FPATH,FNAME]);
[bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr,1);
Scale_mx = eye(2);
Scale_mx(1) = bx(2);
Scale_mx(4) = by(2);
Trans_mx = [bx(1), by(1)]';

parameters.scaler = 1.8;
parameters.small_scaler = 1.8;
parameters.n_points = 9;
parameters.spacing = 175;


% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

if debug_on;
  PsychDebugWindowConfiguration;
end

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if available
eyeTrack_screenNumber = min(screens);
stimulus_screenNumber = max(screens);

% Define black and white
white = WhiteIndex(stimulus_screenNumber);
black = BlackIndex(stimulus_screenNumber);
grey = 128;
inc = white - grey;

% Open an on screen window
eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,grey);
[screenXpixels, screenYpixels] = Screen('WindowSize',eyeTrack_window);
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);
small_rect = [1 1 s1*parameters.small_scaler s2*parameters.small_scaler];

# find unique locations
XXYY = nan * ones(length(trial_records),2);
for i = 1:length(trial_records)
  XXYY(i,1) = trial_records(i).positionX;
  XXYY(i,2) = trial_records(i).positionY;
end

XXYY = unique(XXYY,'rows');

for i = 1:size(XXYY,1);
  small_rects(:,:,i)  = CenterRectOnPoint(small_rect, XXYY(i,1), XXYY(i,2));
end

small_rects = reshape(small_rects,4,9);

calibrating = true;
trans_mtrx = [0;0];
scale_mtrx = eye(2);

dots_colors = [255,0,0;
               0,255,0;
               0,0,255;
               255,255,255;
               0,0,0;
               127,127,127;
               255,255,0;
               0,255,255;
               255,0,255];

while calibrating
  
  Screen('Flip', eyeTrack_window,0,0);
  
  eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage); 
  Screen('DrawTextures', eyeTrack_window, eyeTrack_imageTexture, [], small_rects);
  for t = 1:length(trial_records)
    if strcmp(trial_records(t).trial_error,'no_error')
      Screen('DrawDots',eyeTrack_window,scale_mtrx*trial_records(t).gaze_position+trans_mtrx,3,dots_colors(trial_records(t).pos,:));
    end    
  end

  Screen('Drawtext',eyeTrack_window,'Calibration OK y/n?',2,200)
  Screen('Flip', eyeTrack_window,0,1);
  yn = GetString();
  if yn == 'n'  
    Screen('Drawtext',eyeTrack_window,'Center-to-center spacing is 125px',2,2)
    Screen('Drawtext',eyeTrack_window,'Enter x translation, follow by enter',2,100)
    Screen('Flip', eyeTrack_window,0,1);
    trans_mtrx(1,1) = str2num(GetString());
    Screen('Drawtext',eyeTrack_window,'Enter y translation, follow by enter',2,120)
    Screen('Flip', eyeTrack_window,0,1);
    trans_mtrx(2,1) = str2num(GetString());
    Screen('Drawtext',eyeTrack_window,'Enter x gain, follow by enter',2,140)
    Screen('Flip', eyeTrack_window,0,1);
    scale_mtrx(1,1) = str2num(GetString());
    Screen('Drawtext',eyeTrack_window,'Enter y gain, follow by enter',2,160)
    Screen('Flip', eyeTrack_window,0,0);
    scale_mtrx(2,2) = str2num(GetString());
  elseif yn == 'y'
    calibrating = false;
  end   
end 

Screen('Close',eyeTrack_window);
sca;

# compute new calibration matrices
Scale_mx = Scale_mx.*scale_mtrx;
Trans_mx = Trans_mx + trans_mtrx;

animal = 'Tully-';
saveSTR = [animal,'corrected-calibration-matrices-',date,'.mat'];