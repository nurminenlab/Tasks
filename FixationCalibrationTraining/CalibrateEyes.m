% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;
clear;

% user defined parameters
scaler = 1.2;
stimulus_duration = 2;
reward_size_time = 0.2;

% set-up Datapixx
Datapixx('Open');
Datapixx('StopAllSchedules');
Datapixx('DisableDacAdcLoopback');
Datapixx('RegWrRd');

% allocate space for acquiring eye-tracker voltages
x_center = nan*ones(1,40);
x_off    = nan*ones(1,40);
y_center = nan*ones(1,40);
y_off    = nan*ones(1,40);
n_fixat  = 0;

% Params for acquiring ADC data.
% This streaming demo stores collected ADC data in a 1 second circular buffer within the DATAPixx.
% This circular buffer is then uploaded to a large local matrix at regular intervals.
adcRate = 1e3;                            % Acquire ADC data at 10 kSPS
nAdcSamples = 600;                   % Streaming will use 1 second circular buffer in DATAPixx
adcBuffBaseAddr = 4e6;                      % DATAPixx internal buffer address
minStreamFrames = floor(adcRate / 100);     % Limit streaming reads to 100 times per second


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
[X,Y] = meshgrid(screenXpixels/2-250:250:screenXpixels/2+250,screenYpixels/2-250:250:screenYpixels/2+250);
for i = 1:numel(X);  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
end

% Make the image into a texture
imageTexture = Screen('MakeTexture', window, theImage);

tr_ind = 0;
tr = struct();
while 1  
  [secs, keyCode, deltaSecs] = KbStrokeWait;
  if KbName(keyCode) == '7'        
    pos = 1;
    stim = 'left-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '4'
    pos = 2;
    stim = 'left-center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '1'
    pos = 3;
    stim = 'left-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '8'
    pos = 4;
    stim = 'center-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '5'
    pos = 5;
    stim = 'center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '2'
    pos = 6;
    stim = 'center-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '9'
    pos = 7;
    stim = 'right-top';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '6'
    pos = 8;
    stim = 'right-center';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == '3'
    pos = 9;
    stim = 'right-bottom';
    [rect_center_X,rect_center_Y] = RectCenter(rects(:,:,pos));
  elseif KbName(keyCode) == 'p'
    compute_calibration_matrix(tr,1);
  elseif KbName(keyCode) == 'q'
    break;
  else
    % do nothing
  end
  
  % Draw the image to the screen
  Screen('DrawTexture', window, imageTexture, [], rects(:,:,pos), 0);
  % Start eye tracking and flip the stimulus to the screen  
  Datapixx('SetAdcSchedule', 0, adcRate, nAdcSamples, [0 2]);
  Datapixx('StartAdcSchedule')
  Datapixx('RegWrRd');
  vbl = Screen('Flip', window);
  Screen('FillRect', window, grey, rects(:,:,pos));  
  flipped = 0;
  tic();
  while toc() < stimulus_duration;
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'space';
      Screen('Flip', window);
      flipped = 1;     
      % Wait for the ADC to finish acquiring its scheduled dataset
      while 1
        Datapixx('RegWrRd');   % Update registers for GetAdcStatus
        status = Datapixx('GetAdcStatus');
        if (~status.scheduleRunning)
            break;
        end
      end
    
      % give juice 
      fprintf('juice');
      Datapixx('RegWrRd');
      [adcData, adcTimetags] = Datapixx('ReadAdcBuffer', nAdcSamples);      
      
      # turn pump on
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');
      
      plot(adcTimetags, adcData(1,:)./max(abs(adcData(1,:))),'r-');
      hold on;
      plot(adcTimetags, adcData(2,:)./max(abs(adcData(2,:))),'b-');      
      
      wait_user = 1;
      while wait_user
        [secs, keyCode, deltaSecs] = KbStrokeWait;
        if KbName(keyCode) == 'g'
          tr_ind = tr_ind + 1;
          [X,Y,B] = ginput(2);
          t = find(adcTimetags >= X(1) & adcTimetags <= X(2));      
          tr(tr_ind) = struct('stimulus', stim, "H_voltage", adcData(1,:), "V_voltage", adcData(2,:), ...
                        "H_voltage_median", median(adcData(1,t)), "V_voltage_median", median(adcData(2,t)), ...
                        "rect_center_X", rect_center_X, "rect_center_Y", rect_center_Y);            
          wait_user = 0;
        elseif KbName(keyCode) == 'b'  
          wait_user = 0;
        else
        end                  
      end                
      close;
      break;
    end
  end  
  
  if flipped == 0
    Screen('Flip', window);
  end  
  
end

% Clear the screen
Datapixx('StopAllSchedules');
sca;

%function compute_calibration_matrix(tr)  
%endfunction
