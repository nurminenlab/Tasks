function Texture_task_wIP(debug_on)

% TCP settings to connect with Plotter
t1 = tcp('172.27.85.52', 51007);
t2 = tcp('172.27.85.52', 51008);
  
% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;

# use mouse instead of eye tracker
mouse_track = 1;

% user defined parameters
scaler               = 1.1;
trackWin_factor      = 1.3;
wait_fixation        = 4;
rewardConsume_period = 2;
max_fixation_time    = 4;
ms                   = 10;
min_target_time      = 0.5;
max_trs              = 10000;
response_wait_min    = 0.125;
response_wait_max    = 0.225;
gaze_move_time       = 10;

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
parameters.contrast_C = 1;

thetas = deg2rad([-45,45,135,225]);
R      = 200;

fix_point_Window_size = 100;
trackMarkerColor = [255,0,0];

if mouse_track
  XY = ones(2,1)*nan;
end

Datapixx('Open')
adcRate = 1e3;
baseBuffAddr = 0;
minStreamFrames = 15;

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
stimulus_window = Screen('OpenWindow',stimulus_screenNumber,grey);
eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,grey);
windowPointer = stimulus_window;

% Get the size of the on screen window, these are the same for both screens
[screenXpixels, screenYpixels] = Screen('WindowSize',stimulus_window);

% Load marmoset face
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);

% scale image rectangle
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

% Make the image into a texture
stimulus_imageTexture = Screen('MakeTexture', stimulus_window, theImage);
eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);

fix_point_rect = [1 1 10 10];
fix_point_rect = CenterRectOnPoint(fix_point_rect, screenXpixels/2, screenYpixels/2);

fix_point_Window = fix_point_Window_size;
fix_point_Windowrect = [1 1 fix_point_Window fix_point_Window];
fix_point_Windowrect = CenterRectOnPoint(fix_point_Windowrect,screenXpixels/2, screenYpixels/2);
fix_point_Window = fix_point_Window/2;

target_diameter = parameters.plateauPixels_C + parameters.edgePixels_C;
target_Windowrect = [1 1 target_diameter target_diameter];
target_Windowrect = CenterRectOnPoint(target_Windowrect,screenXpixels/2, screenYpixels/2);

txt_rect = [1 1 gridSize gridSize];
txt_rect = CenterRectOnPoint(txt_rect, screenXpixels/2, screenYpixels/2);
txt_rects  = nan * ones(4,length(thetas)*length(R));

[Xoff, Yoff] = pol2cart(thetas,R);
for i = 1:length(Xoff)
  txt_rects(:,i) = round(CenterRectOnPoint(txt_rect, screenXpixels/2+Xoff(i), screenYpixels/2+Yoff(i)));
end

[grText_CS,grText_S,grText_CS_iTrack,grText_S_iTrack] = generate_filteredNoise_CS(gridSize,stimulus_window,eyeTrack_window,parameters);

tr_ind = 0;
tr = struct();
is_running = 1;
conditions = ['5'];
pos = 5;
stim = 'center';

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

  tr_ind = tr_ind + 1
  
  waiting_for_response = 0;
  tracking_reward = 0;
  on_target = 0;
  
  # target position
  target_pos = randi(4);
  min_target_time = response_wait_min + (response_wait_max - response_wait_min)*rand();
  
  if tr_ind > 1
    Screen('Flip', eyeTrack_window,0);
  end   
  
  % Grey screen
  Screen('FillRect', stimulus_window, grey, rects(:,:,pos));    
  Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));  
  greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);  
  greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window);  
  
  target_Windowrect = round(CenterRectOnPoint(target_Windowrect, screenXpixels/2+Xoff(target_pos), screenYpixels/2+Yoff(target_pos)));
  target_pos_X = screenXpixels/2+Xoff(target_pos);  
  target_pos_Y = screenYpixels/2+Yoff(target_pos); 
  
  not_target_Xoff = Xoff;
  not_target_Xoff(target_pos) = [];  
  not_target_pos_X = screenXpixels/2+not_target_Xoff;
  
  not_target_Yoff = Yoff;
  not_target_Yoff(target_pos) = [];
  not_target_pos_Y = screenYpixels/2+not_target_Yoff;
  
  grText = repmat(grText_S,1,4);
  grText(:,target_pos) = grText_CS;
  grText_iTrack = repmat(grText_S_iTrack,1,4);
  grText_iTrack(:,target_pos) = grText_CS_iTrack;
  
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
      
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
        tcp_write(t2,"q");
        is_running = 0;      
        close all;
        clear t1;
        clear t2;
        sca;
        Datapixx('StopAllSchedules');
        Datapixx('Close');
        break;
    end    
  end  %
    
  % Draw monkey face
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));    
  % Draw fixation window
  Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
  stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period);
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
      
      Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));           
      Screen('FillOval', eyeTrack_window,[0 0 255],fix_point_rect);       
      Screen('DrawTextures', eyeTrack_window, grText_iTrack,[],txt_rects);
      Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3);
      
      Screen('FillRect', stimulus_window, grey, rects(:,:,pos));           
      Screen('FillOval', stimulus_window,[0 0 255],fix_point_rect);
      Screen('DrawTextures', stimulus_window, grText,[],txt_rects);      

      # write(s, "stimulus information")
      % sending the stimlus condition 
      disp(target_pos);
      tcp_write(t1,num2str(target_pos));
          
      # gaze on monkey face, draw targets on screen, start recording
      tcp_write(t2,"1");
      
      Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);     
      Screen('Flip', eyeTrack_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
      break;
    end
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    reaction_time = nan;
    trial_error = 'no_fixation';
    on_target_time = nan;
    
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
        tcp_write(t2,"q");
        is_running = 0;      
        close all;
        clear t1;
        clear t2;
        sca;
        Datapixx('StopAllSchedules');
        Datapixx('Close');
        break;
    end  
    
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
      trial_error = 'broke_fixation';
      # stop recording
      break;      
    end           
   
    # the animals needs to hold fixation for "min_target_time" before responding
    if toc(eyeTrack_clock) > min_target_time
      on_target = 0;
      waiting_for_response = 1;
      on_target_time = toc(eyeTrack_clock);
      trial_error = 'no_error';
      Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));
      Screen('FrameOval',eyeTrack_window,[0 0 255],target_Windowrect,3,3);
      Screen('FillRect', stimulus_window, grey, rects(:,:,pos));
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      response_time_clock = tic();      
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);   
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
        tcp_write(t2,"q");
        is_running = 0;      
        close all;
        clear t1;
        clear t2;
        sca;
        Datapixx('StopAllSchedules');
        Datapixx('Close');
        break;
    end        
  end  %
  
  gaze_moved = 0;
  while waiting_for_response;     
    
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
   
    # if eyes left fixation window the monkey needs to go to the target in an x amount of time
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow
      gaze_move_clock = tic();
      gaze_moved      = 1;
    end
    
    if gaze_moved && (toc(gaze_move_clock) > gaze_move_time)
      # stop recording
      tcp_write(t2,"0");
      break;
    end
    
    # target hit
    if sqrt((XY(1) - target_pos_X)^2 + (XY(2) - target_pos_Y)^2) < 150
      Screen('FillRect', eyeTrack_window, grey, txt_rects);
      Screen('FillRect', stimulus_window, grey, txt_rects);
      % Draw monkey face
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], txt_rects(:,target_pos));  
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], txt_rects(:,target_pos));
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      reward_time_clock = tic();
      tracking_reward = 1;
      # stop recording
      tcp_write(t2,"0");
      break;
    end
    
##    if sqrt((XY(1) - not_target_pos_X(1))^2 + (XY(2) - not_target_pos_Y(1))^2) < 150
##      waiting_for_response = 0;
##      break;
##    end
##    
##    if sqrt((XY(1) - not_target_pos_X(2))^2 + (XY(2) - not_target_pos_Y(2))^2) < 150
##      waiting_for_response = 0;
##      break;
##    end
##    
##    if sqrt((XY(1) - not_target_pos_X(3))^2 + (XY(2) - not_target_pos_Y(3))^2) < 150
##      waiting_for_response = 0;
##      break;
##    end
##        
##    # max length of trial
##    if toc(response_time_clock) > 3      
##      waiting_for_response = 0;      
##      break;
##    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);    
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
        tcp_write(t2,"q");
        is_running = 0;      
        close all;
        clear t1;
        clear t2;
        sca;
        Datapixx('StopAllSchedules');
        Datapixx('Close');
        break;
    end    
  end  %

  while tracking_reward;
    
    if ~mouse_track
      Datapixx('RegWrRd');if keyIsDown && KbName(keyCode) == 'q';
      tcp_write(t2,"q");
      is_running = 0;      
      close all;
      clear t1;
      clear t2;
      sca;
      Datapixx('StopAllSchedules');
      Datapixx('Close');
      break;
    end
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
   
    if sqrt((XY(1) - target_pos_X)^2 + (XY(2) - target_pos_Y)^2) > 100
      Screen('FillRect', eyeTrack_window, grey);
      Screen('FillRect', stimulus_window, grey);  
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      tracking_reward = 0;
      
      reward_size_time = 0.25*sqrt((toc(reward_time_clock)));
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();      
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');      
      break;
    end
    
    # max length of trial
    if toc(reward_time_clock) > 3
      reward_size_time = 0.25*sqrt((toc(reward_time_clock)));
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();      
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');
      tracking_reward = 0;      
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);     
  end  %  
  
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      clear t1;
      clear t2;
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
clear t1;
clear t2;