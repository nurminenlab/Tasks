% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;
clear;

# use mouse instead of eye tracker
mouse_track = 1;

on_target = 0;

% user defined parameters
scaler = 1.1;
wait_fixation = 4;
rewardConsume_period = 2;
max_fixation_time = 4;
ms = 10;
min_target_time = 0.15;
max_trs = 10000;

if mouse_track
  XY = ones(2,1)*nan;
end

Datapixx('Open')
adcRate = 1e3;
baseBuffAddr = 0;
%nAdcSamples = floor((wait_fixation + rewardConsume_period + max_fixation_time) * adcRate);
minStreamFrames = 15;

trackMarkerColor = [255,0,0];

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

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

#PsychDebugWindowConfiguration;
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
#[stimulus_window, stimulus_windowRect] = PsychImaging('OpenWindow', stimulus_screenNumber, grey);
#[eyeTrack_window, eyeTrack_windowRect] = PsychImaging('OpenWindow', eyeTrack_screenNumber, grey);
stimulus_window = Screen('OpenWindow',stimulus_screenNumber,128);
eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,128);

gridSize = 128;
orientations = [0,90];
pixelsPerPeriod = 33;
plateauCycles = 3;
edgeCycles = 0.25;
contrast = 0.3
windowPointer = stimulus_window;

% Get the size of the on screen window, these are the same for both screens
[screenXpixels, screenYpixels] = Screen('WindowSize',stimulus_window);

% Set up alpha-blending for smooth (anti-aliased) lines
#Screen('BlendFunction', stimulus_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
#Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Load marmoset face
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);

% scale image rectangle
trackWin_factor = 1.3;
rect = [0 0 s1*scaler s2*scaler];
eyePos_rect = [0 0 5 5];
trackWindow_rect = [0 0 s1*scaler*trackWin_factor s2*scaler*trackWin_factor];
trackWindow = s1*scaler*trackWin_factor/2;

idx = 0;
XBC = [screenXpixels/2-125, screenXpixels/2, screenXpixels/2 + 125];
YBC = [screenYpixels/2-125, screenYpixels/2, screenYpixels/2 + 125];
[X,Y] = meshgrid(XBC,YBC);

for i = 1:numel(X);  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
  trackWindow_rect(:,:,i) = CenterRectOnPoint(trackWindow_rect, X(i), Y(i));
end

grating_rect = [0 0 gridSize gridSize];
grating_rect = CenterRectOnPoint(grating_rect,screenXpixels/2-135,screenYpixels/2+135);

% Make the image into a texture
stimulus_imageTexture = Screen('MakeTexture', stimulus_window, theImage);
eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);

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

txt_rect = [1 1 gridSize gridSize];
txt_rect = CenterRectOnPoint(txt_rect, screenXpixels/2, screenYpixels/2);
txt_rects  = nan * ones(4,4);
thetas = deg2rad([-45,45,135,225]);
R = 200;
[Xoff, Yoff] = pol2cart(thetas,R);
for i = 1:length(Xoff)
  txt_rects(:,i) = round(CenterRectOnPoint(txt_rect, screenXpixels/2+Xoff(i), screenYpixels/2+Yoff(i)));
end

[grText_CS,grText_S,grText_CS_iTrack,grText_S_iTrack] = generate_filteredNoise_CS(gridSize,stimulus_window,eyeTrack_window,parameters);
grText        = [repmat(grText_S,1,3),grText_CS];
grText_iTrack = [repmat(grText_S_iTrack,1,3),grText_CS_iTrack];

tr_ind = 0;
tr = struct();
#conditions = ['1','2','3','4','5','6','7','8','9'];
#conditions = ['4','5','6','8','2'];
conditions = ['5'];
is_running = 1;

KbStrokeWait();
if ~mouse_track
  Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
  Datapixx('StartAdcSchedule');
  Datapixx('RegWrRd');
  # to force correct baseBuffAddr
  XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
  Datapixx('RegWrRd');
end

while is_running
  tr_ind = tr_ind + 1;
  %eyePos_rect = CenterRectOnPoint(eyePos_rect, 960, 540-150);

  if tr_ind > 1
    Screen('Flip', eyeTrack_window,0);
  end 
  
  % randomize conditon
  c = conditions(randi(length(conditions)));  
  if c == '7'        
    pos = 1;
    stim = 'left-top';    
  elseif c == '4'
    pos = 2;
    stim = 'left-center';    
  elseif c == '1'
    pos = 3;
    stim = 'left-bottom';    
  elseif c == '8'
    pos = 4;
    stim = 'center-top';    
  elseif c == '5'
    pos = 5;
    stim = 'center';    
  elseif c == '2'
    pos = 6;
    stim = 'center-bottom';    
  elseif c == '9'
    pos = 7;
    stim = 'right-top';    
  elseif c == '6'
    pos = 8;
    stim = 'right-center';    
  elseif c == '3'
    pos = 9;
    stim = 'right-bottom';    
  else
    % do nothing
  end
  
  % Grey screen
  Screen('FillRect', stimulus_window, grey, rects(:,:,pos));  
  Screen('FillRect', stimulus_window, grey, grating_rect);
  Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));  
  greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);
  #greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window,0,1);
  greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window);
  
  tic();
  while toc() < rewardConsume_period;
      [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();                    
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
        [x,y,buttons,focus] = GetMouse(eyeTrack_window);
        XY(1) = x;
        XY(2) = y;
      end      
      
      eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
      Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);      
      Screen('Flip', eyeTrack_window,0,1);    
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;
      Datapixx('StopAllSchedules');
      Datapixx('Close');      
      close all;
      sca;
      break;
    end
  end  %
    
  % Draw monkey face
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));    
  % Draw fixation window
  Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
  stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
  stimulusScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window, greyScreen_eyeTrack_vbl + rewardConsume_period,1);
  
  wait_fixation_clock = tic();
  while toc(wait_fixation_clock) < wait_fixation;    
    
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
      [x,y,buttons,focus] = GetMouse(eyeTrack_window);
      XY(1) = x;
      XY(2) = y;
    end           
    
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) < trackWindow
      on_target = 1;
      reaction_time = toc(wait_fixation_clock);
      eyeTrack_clock = tic();      
      Screen('DrawTextures', stimulus_window, grText,[],txt_rects);      
      Screen('DrawTextures', eyeTrack_window, grText_iTrack,[],txt_rects);
      Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period);     
      Screen('Flip', eyeTrack_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
      break;
    end
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    reaction_time = nan;
    trial_error = 'no_fixation';
    on_target_time = nan;
    
  end  %

  while on_target;     
    
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
      [x,y,buttons,focus] = GetMouse(eyeTrack_window);
      XY(1) = x;
      XY(2) = y;
    end      
    
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow
      on_target = 0;
      on_target_time = toc(eyeTrack_clock);
      if on_target_time < min_target_time
        trial_error = 'broke_fixation';               
      else
        trial_error = 'no_error';        
        # give juice
        # turn pump on
        Datapixx('SetDoutValues', 1);
        Datapixx('RegWrRd');
        a = tic();
        reward_size_time = 0.25*sqrt((on_target_time));
        while toc(a) < reward_size_time
          # pump juice
        end
        # turn pump off
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');        
      end      
      break;
    end 
   
    # max length of trial
    if toc(eyeTrack_clock) > max_fixation_time
      on_target = 0;
      on_target_time = toc(eyeTrack_clock);
      trial_error = 'no_error';
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();
      reward_size_time = 0.3*sqrt((on_target_time));
      while toc(a) < reward_size_time
         # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');      
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      break;
    end    
  end  %
  
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      break;
  end
 
  % record trial data
  trial_records(tr_ind).stimulus = stimulus_image;
  trial_records(tr_ind).positionX = X(pos);
  trial_records(tr_ind).positionY = Y(pos);
  trial_records(tr_ind).reaction_time  = reaction_time;
  trial_records(tr_ind).on_target_time = on_target_time; 
  trial_records(tr_ind).trial_error = trial_error; 
  
  if tr_ind >= max_trs
    is_running = 0;
    sca;
  end
  
end


