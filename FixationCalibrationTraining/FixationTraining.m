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

# load fixation calibration
[FNAME, FPATH, FLTIDX] = uigetfile();
load([FPATH,FNAME]);
[bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr);
Scale_mx = eye(2);
Scale_mx(1) = bx(2);
Scale_mx(4) = by(2);
Trans_mx = [bx(1), by(1)]';

% set-up Datapixx
Datapixx('Open');
adcRate = 1e3;
nAdcSamples = floor((stimulus_duration + isi_duration) * adcRate);
minStreamFrames = 15;

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

% Get the size of the on screen window, these are the same for both screens
[screenXpixels, screenYpixels] = Screen('WindowSize', stimulus_window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(stimulus_windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', stimulus_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Load marmoset face
theImage = imread('face8.jpg');
[s1, s2, s3] = size(theImage);

% scale image rectangle
rect = [0 0 s1*scaler s2*scaler];

idx = 0;
[X,Y] = meshgrid(screenXpixels/2-250:250:screenXpixels/2+250,screenYpixels/2-250:250:screenYpixels/2+250);
for i = 1:numel(X);  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
end

% Make the image into a texture
stimulus_imageTexture = Screen('MakeTexture', stimulus_window, theImage);
eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);

tr_ind = 0;
tr = struct();
#conditions = ['1','2','3','4','5','6','7','8','9'];
conditions = ['4','5','6'];
is_running = 1;

KbStrokeWait();

while is_running
  tr_ind = tr_ind + 1;
  
  c = conditions(randi(length(conditions)));  
  if c == '7'        
    pos = 1;
    stim = 'left-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '4'
    pos = 2;
    stim = 'left-center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '1'
    pos = 3;
    stim = 'left-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '8'
    pos = 4;
    stim = 'center-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '5'
    pos = 5;
    stim = 'center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '2'
    pos = 6;
    stim = 'center-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '9'
    pos = 7;
    stim = 'right-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '6'
    pos = 8;
    stim = 'right-center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif c == '3'
    pos = 9;
    stim = 'right-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));  
  else
    % do nothing
  end
  
  % Grey screen
  Screen('FillRect', stimulus_window, grey, rects(:,:,pos));  
  Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));
  stimulus_vbl = Screen('Flip', stimulus_window);
  eyeTrack_vbl = Screen('Flip', eyeTrack_window);
  
  % Draw monkey face
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos), 0);
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos), 0);  
  Screen('Flip', stimulus_window, stimulus_vbl + isi_duration,0);
  Screen('Flip', eyeTrack_window, eyeTrack_vbl + isi_duration,1);
  
  Datapixx('SetAdcSchedule', 0, adcRate, nAdcSamples, [0 2]);
  Datapixx('StartAdcSchedule')  
  tic();
  while toc() < isi_duration;
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    % Start eye tracking and flip the stimulus to the screen        
    Datapixx('RegWrRd');
    status = Datapixx('GetAdcStatus');
    nReadFrames = status.newBufferFrames;   
    if (nReadFrames < minStreamFrames)
      continue;
    else
      dada = Datapixx('ReadAdcBuffer', nReadFrames, -1);   
      dada = Scale_mx*dada + Trans_mx;
      Screen('DrawLines', eyeTrack_window, dada);    
      Screen('Flip', eyeTrack_window,0,1);
    end        
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      Datapixx('StopAllSchedules');
      Datapixx('Close');
      sca;
      break;
    end
  end
  
  tic();
  while toc() < stimulus_duration;
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    % Start eye tracking and flip the stimulus to the screen        
    Datapixx('RegWrRd');
    status = Datapixx('GetAdcStatus');
    nReadFrames = status.newBufferFrames;   
    if (nReadFrames < minStreamFrames)
      continue;
    else
      dada = Datapixx('ReadAdcBuffer', nReadFrames, -1);   
      dada = Scale_mx*dada + Trans_mx;
      Screen('DrawLines', eyeTrack_window, dada);    
      Screen('Flip', eyeTrack_window,0,1);
    end        
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      Datapixx('StopAllSchedules');
      Datapixx('Close');
      sca;
      break;
    end
  end  
  
  Screen('Flip', eyeTrack_window,0);
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      Datapixx('StopAllSchedules');
      Datapixx('Close');
      sca;
      break;
  end
  
end


