% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;
clear;

% user defined parameters
scaler = 3.5;
stimulus_duration = 10;
isi_duration = 10;

reward_size_time = 0.2;
inter_juice_time = 1-reward_size_time;

% set-up Datapixx
##Datapixx('Open');
##Datapixx('StopAllSchedules');
##Datapixx('DisableDacAdcLoopback');
##Datapixx('RegWrRd');

% allocate space for acquiring eye-tracker voltages
x_center = nan*ones(1,40);
x_off    = nan*ones(1,40);
y_center = nan*ones(1,40);
y_off    = nan*ones(1,40);
n_fixat  = 0;

Datapixx('Open');
Datapixx('StopAllSchedules');
Datapixx('RegWrRd'); 

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if available
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Load marmoset face
theImage = imread('face8.jpg');
[s1, s2, s3] = size(theImage);

% scale image rectangle
rect = [0 0 s1*scaler s2*scaler];

idx = 0;
%[X,Y] = meshgrid(screenXpixels/2-250:250:screenXpixels/2+250,screenYpixels/2-250:250:screenYpixels/2+250);
[X,Y] = meshgrid(screenXpixels/2-250:250:screenXpixels/2+250,screenYpixels/2+250);
npos = numel(X);
for i = 1:npos;  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
end

% Make the image into a texture
imageTexture = Screen('MakeTexture', window, theImage);

tr_ind = 0;
tr = struct();
first_time = true;
max_trials = 500;
n_drops = 9;

while tr_ind < max_trials;
  tr_ind = tr_ind + 1;
  FlushEvents('keyDown');
  #[secs, keyCode, deltaSecs] = KbStrokeWait;
  pos = randi(npos);  
  [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  
  
  % Draw the image to the screen
  Screen('DrawTexture', window, imageTexture, [], rects(:,:,pos), 0);
  % Start eye tracking and flip the stimulus to the screen  
##  Datapixx('SetAdcSchedule', 0, adcRate, nAdcSamples, [0 2]);
##  Datapixx('StartAdcSchedule')
##  Datapixx('RegWrRd');
  if first_time
    vbl = Screen('Flip', window);
    first_time = false;
  else
    vbl = Screen('Flip', window, vbl+isi_duration);
  end
  
  # we give n juice drops
  for i = 1:n_drops
    # turn pump on
    Datapixx('SetDoutValues', 1);
    Datapixx('RegWrRd');
    a = tic();
    while toc(a) < reward_size_time
      # do nothing
    end
    # turn pump off
    Datapixx('SetDoutValues', 0);
    Datapixx('RegWrRd');    
    a = tic();
    while toc(a) < inter_juice_time
      # do nothing
    end
  end
  
  
  Screen('FillRect', window, grey, rects(:,:,pos));  
  vbl = Screen('Flip', window, vbl+stimulus_duration);  
  
end

% Clear the screen
#Datapixx('StopAllSchedules');
sca;

%function compute_calibration_matrix(tr)  
%endfunction
