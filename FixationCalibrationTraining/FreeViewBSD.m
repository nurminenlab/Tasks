% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;
clear;

debug_on = 1;
mouse_track = 0;

adcRate = 1e3;
baseBuffAddr = 0;
%nAdcSamples = floor((wait_fixation + rewardConsume_period + max_fixation_time) * adcRate);
minStreamFrames = 15;

if debug_on;
  PsychDebugWindowConfiguration;
end

imdir = '/home/vpixx/Tasks/BSDS300/images/all';
fls = dir(imdir);
fls([1,2]) = [];

% user defined parameters
scaler = 2;
stimulus_duration = 10;
isi_duration = 10;
trackMarkerColor = [255,0,0];

reward_size_time = 0.2;
inter_juice_time = 1-reward_size_time;

% allocate space for acquiring eye-tracker voltages
x_center = nan*ones(1,40);
x_off    = nan*ones(1,40);
y_center = nan*ones(1,40);
y_off    = nan*ones(1,40);
n_fixat  = 0;

if mouse_track
  XY = ones(2,1)*nan;
end

Datapixx('Open');
Datapixx('StopAllSchedules');
Datapixx('RegWrRd'); 

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if available
consoleScreen   = min(screens);
projectorScreen = max(screens);

% Define black and white
white = WhiteIndex(projectorScreen);
black = BlackIndex(projectorScreen);
grey = white / 2;
inc = white - grey;

% Open an on screen window
[projectorWindow, windowRect_projector] = PsychImaging('OpenWindow', projectorScreen, grey);
[consoleWindow, windorRect_console]     = PsychImaging('OpenWindow', consoleScreen, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', projectorWindow);

% Query the frame duration
ifi = Screen('GetFlipInterval', projectorWindow);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect_projector);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', projectorWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('BlendFunction', consoleWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Load marmoset face
#fl_inds  = randperm(length(fls));
fl_inds  = [1,2];
n_images = 2;
imageTexture = nan * ones(1,n_images);
imageTexture_console = nan * ones(1,n_images);

for i = 1:n_images
  % Make the image into a texture
  theImage = imread(fullfile(fls(fl_inds(i)).folder,fls(fl_inds(i)).name));
  imageTexture(i) = Screen('MakeTexture', projectorWindow, theImage);
  imageTexture_console(i) = Screen('MakeTexture', consoleWindow, theImage);
end

[s1, s2, s3] = size(theImage);
scaler = ceil(screenXpixels/s2);
% scale image rectangle
rect = [0 0 s2*scaler s1*scaler];
rect = CenterRectOnPoint(rect,screenXpixels/2,screenYpixels/2);
eyePos_rect = [0 0 5 5];

tr_ind = 0;
tr = struct();
first_time = true;
max_trials = 500;
n_drops = 20;

if ~mouse_track
  # load fixation calibrationscreenXpixels/2
  [FNAME, FPATH, FLTIDX] = uigetfile();
  load([FPATH,FNAME]);
  [bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr,1);
  Scale_mx = eye(2);
  Scale_mx(1) = bx(2);
  Scale_mx(4) = by(2);
  Trans_mx = [bx(1), by(1)]';
end

KbStrokeWait();
if ~mouse_track
  Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
  Datapixx('StartAdcSchedule');
  Datapixx('RegWrRd');
  # to force correct baseBuffAddr
  XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
  Datapixx('RegWrRd');
end

while tr_ind < max_trials;
  tr_ind = tr_ind + 1;  
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck(); 
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;
      Datapixx('StopAllSchedules');
      Datapixx('Close');      
      close all;
      sca;
      #save(saveSTR,'trial_records');
      break;
  end
  
  i_ind = rem(tr_ind,2) + 1;  
  % Draw the image to the screen
  Screen('DrawTexture', projectorWindow, imageTexture(i_ind), [], rect);
  Screen('DrawTexture', consoleWindow, imageTexture_console(i_ind), [], rect);
  
  % Start eye tracking and flip the stimulus to the screen  
  if ~mouse_track
    Datapixx('RegWrRd');
    status = Datapixx('GetAdcStatus');
    nReadFrames = status.newBufferFrames;   
    if (nReadFrames < minStreamFrames)
      continue;
    else        
      XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
      XY = Scale_mx*XY + Trans_mx;
    end
  else
    [x,y,buttons,focus] = GetMouse(consoleWindow);
    XY(1) = x;
    XY(2) = y;
  end
  
  eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
  Screen('FillOval', consoleWindow, trackMarkerColor, eyePos_rect);
 
  Screen('Flip', projectorWindow);
  Screen('Flip', consoleWindow,0,1);
  
  # we give n juice drops
  for i = 1:n_drops
    # turn pump on
    Datapixx('SetDoutValues', 1);
    Datapixx('RegWrRd');
    a = tic();
    while toc(a) < reward_size_time
      [keyIsDown, secs, keyCode, deltaSecs] = KbCheck(); 
      if keyIsDown && KbName(keyCode) == 'q';
        is_running = 0;
        Datapixx('StopAllSchedules');
        Datapixx('Close');      
        close all;
        sca;
        #save(saveSTR,'trial_records');
        break;
      end
      
      if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else        
          XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
          XY = Scale_mx*XY + Trans_mx;
        end
      else
        [x,y,buttons,focus] = GetMouse(consoleWindow);
        XY(1) = x;
        XY(2) = y;
      end  
      eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
      Screen('FillOval', consoleWindow, trackMarkerColor, eyePos_rect);
      Screen('Flip', consoleWindow,0,1);
    end
    
    # turn pump off
    Datapixx('SetDoutValues', 0);
    Datapixx('RegWrRd');    
    a = tic();
    while toc(a) < inter_juice_time      
      [keyIsDown, secs, keyCode, deltaSecs] = KbCheck(); 
      if keyIsDown && KbName(keyCode) == 'q';
        is_running = 0;
        Datapixx('StopAllSchedules');
        Datapixx('Close');      
        close all;
        sca;
        #save(saveSTR,'trial_records');
        break;
      end
      
      # do nothing
      if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else        
          XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
          XY = Scale_mx*XY + Trans_mx;
        end
      else
        [x,y,buttons,focus] = GetMouse(consoleWindow);
        XY(1) = x;
        XY(2) = y;
      end  
      eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
      Screen('FillOval', consoleWindow, trackMarkerColor, eyePos_rect);
      Screen('Flip', consoleWindow,0,1);
    end
      
  end  
  
end

sca;

