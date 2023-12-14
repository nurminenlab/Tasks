% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;
clear;

on_target = 0;

imdir = '/home/vpixx/Tasks/MonkeyFaces';
fls = dir(imdir);
fls([1,2,3]) = [];

% user defined parameters
scaler = 1.2;
wait_fixation = 1;
rewardConsume_period = 0.4;
max_fixation_time = 20;
ms = 10;
min_target_time = 0.15;

% set-up Datapixx
Datapixx('Open')
adcRate = 1e3;
baseBuffAddr = 0;
%nAdcSamples = floor((wait_fixation + rewardConsume_period + max_fixation_time) * adcRate);
minStreamFrames = 15;

trackMarkerColor = [255,0,0];
# load fixation calibration
[FNAME, FPATH, FLTIDX] = uigetfile();
load([FPATH,FNAME]);
[bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr,1);
Scale_mx = eye(2);
Scale_mx(1) = bx(2);
Scale_mx(4) = by(2);
Trans_mx = [bx(1), by(1)]';

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
grey = white / 2;
inc = white - grey;

% Open an on screen window
[stimulus_window, stimulus_windowRect] = PsychImaging('OpenWindow', stimulus_screenNumber, grey);
[eyeTrack_window, eyeTrack_windowRect] = PsychImaging('OpenWindow', eyeTrack_screenNumber, grey);

% Get the size of the on screen window, these are the same for both screens
[screenXpixels, screenYpixels] = Screen('WindowSize', stimulus_window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(stimulus_windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', stimulus_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Load marmoset face
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);

% scale image rectangle
trackWin_factor = 2;
rect = [0 0 s1*scaler s2*scaler];
eyePos_rect = [0 0 5 5];
trackWindow_rect = [0 0 s1*scaler*trackWin_factor s1*scaler*trackWin_factor];
trackWindow = s1*scaler*trackWin_factor/2;

idx = 0;
XBC = [screenXpixels/2-125, screenXpixels/2, screenXpixels/2 + 125];
YBC = [screenYpixels/2-125, screenYpixels/2, screenYpixels/2 + 125];
[X,Y] = meshgrid(XBC,YBC);

for i = 1:numel(X);  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
  trackWindow_rect(:,:,i) = CenterRectOnPoint(trackWindow_rect, X(i), Y(i));
end

% Make the image into a texture
num_conds = length(fls);
for i = 1:num_conds
  theImage = imread(fullfile(fls(i).folder,fls(i).name));  
  stimulus_imageTexture(i) = Screen('MakeTexture', stimulus_window, theImage);
  eyeTrack_imageTexture(i) = Screen('MakeTexture', eyeTrack_window, theImage);
end

tr_ind = 0;
tr = struct();
conditions = ['4','5','6','2','8'];
conditions = ['5'];
is_running = 1;

KbStrokeWait();
Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
Datapixx('StartAdcSchedule');
Datapixx('RegWrRd');

# to force correct baseBuffAddr
XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
Datapixx('RegWrRd');
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
  Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));
  greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);
  greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window,0,1);   
  
  tic();
  while toc() < rewardConsume_period;
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();                    
      
      Datapixx('RegWrRd');
      status = Datapixx('GetAdcStatus');
      nReadFrames = status.newBufferFrames;   
      if (nReadFrames < minStreamFrames)
        continue;
      else        
        XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
        XY = Scale_mx*XY + Trans_mx;
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
  
  hold on 
  % Draw monkey face
  im_num = randi(num_conds);
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture(im_num), [], rects(:,:,pos), 0);
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture(im_num), [], rects(:,:,pos), 0);  
  % Draw fixation window
  Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
  stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period);
  stimulusScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window, greyScreen_eyeTrack_vbl + rewardConsume_period,1);
  
  wait_fixation_clock = tic();
  while toc(wait_fixation_clock) < wait_fixation;    
         
    Datapixx('RegWrRd');
    status = Datapixx('GetAdcStatus');%PsychDebugWindowConfiguration;
    nReadFrames = status.newBufferFrames;   
    if (nReadFrames < minStreamFrames)
      continue;
    else      
      XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
      XY = Scale_mx*XY + Trans_mx;
    end   
    
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) < trackWindow
      on_target = 1;
      reaction_time = toc(wait_fixation_clock);
      eyeTrack_clock = tic();      
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
    
    Datapixx('RegWrRd');
    status = Datapixx('GetAdcStatus');
    nReadFrames = status.newBufferFrames;   
    if (nReadFrames < minStreamFrames)
      continue;
    else
      XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
      XY = Scale_mx*XY + Trans_mx;
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
        reward_size_time = 0.8*sqrt((on_target_time));
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
      reward_size_time = 0.4*sqrt((on_target_time));
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
  
end


