function trial_records = RF_map(debug_on)
% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;

PsychDefaultSetup(2);

# use mouse instead of eye tracker
mouse_track = 0;

animal = 'Sansa-';
saveSTR = [animal,'RF-map-trial_records-',date,'.mat'];
save_append = 0;
while exist(saveSTR,'file') == 2
  save_append = save_append + 1;
  saveSTR = [animal,'RF-map-trial_records-',date,'_v',num2str(save_append),'.mat'];
end  

% user defined parameters
if strcmp(animal,'Apollo-')
  scaler               = 1.8;
  small_scaler         = 1.8;
  trackWin_factor      = 2.0;
  wait_fixation        = 0.75;
  rewardConsume_period = 2;
  ms                   = 10;
  min_target_time      = 0.5;
  max_trs              = 10000;
  response_wait_min    = 0.125;
  response_wait_max    = 1;
  gaze_move_time       = 1;
  max_fixation_time    = 1;
  min_fixation_time    = 0.1;
  reward_scaler = 0.6;
  FR = 120;
  gridSize = 256;
  fix_point_Window_size = 100;
  trackMarkerColor = [255,0,0];
  gaze_position = nan*ones(2,FR*ceil((wait_fixation+max_fixation_time)));
  
elseif strcmp(animal,'Sansa-')
  scaler               = 0.7;
  trackWin_factor      = 5.5;
  wait_fixation        = 0.75;
  rewardConsume_period = 2;
  ms                   = 10;
  min_target_time      = 0.5;
  max_trs              = 10000;
  response_wait_min    = 0.125;
  response_wait_max    = 1;
  gaze_move_time       = 1;
  max_fixation_time    = 4;
  min_fixation_time    = 0.1;
  reward_scaler = 0.9;
  FR = 120;
  gridSize = 256;
  fix_point_Window_size = 100;
  trackMarkerColor = [255,0,0];
  gaze_position = nan*ones(2,FR*ceil((wait_fixation+max_fixation_time)));
else
  fprintf('Strange animal \n');
end

  
if mouse_track
  XY = ones(2,1)*nan;
end

Datapixx('Open')
Datapixx('StopAllSchedules')
Datapixx('DisableDacAdcLoopback');
Datapixx('RegWrRd');
adcRate = 1e3;
baseBuffAddr = 5e6;
minStreamFrames = 15;

if ~mouse_track
  # load fixation calibrationscreenXpixels/2
  [FNAME, FPATH, FLTIDX] = uigetfile();
  load([FPATH,FNAME]);
  [bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr,1);
  Scale_mx = eye(2);
  Scale_mx(1) = bx(2);
  Scale_mx(4) = by(2);
  Trans_mx = [bx(1), by(1)-65]';
end


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
stimulus_image = 'face10.jpg';
theImage = imread(stimulus_image);

[s1, s2, s3] = size(theImage);

% scale image rectangle
rect = [0 0 s1*scaler s2*scaler];
#small_rect = [0 0 s1*small_scaler s2*small_scaler];

eyePos_rect = [0 0 5 5];
trackWindow_rect = [0 0 s1*scaler*trackWin_factor s2*scaler*trackWin_factor];
trackWindow = s1*scaler*trackWin_factor/2;

idx = 0;
XBC = [screenXpixels/2-125, screenXpixels/2, screenXpixels/2 + 125];
YBC = [screenYpixels/2-125, screenYpixels/2, screenYpixels/2 + 125];
[X,Y] = meshgrid(XBC,YBC);

for i = 1:numel(X);  
  rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
  #small_rects(:,:,i)  = CenterRectOnPoint(rect, X(i), Y(i));
  trackWindow_rect(:,:,i) = CenterRectOnPoint(trackWindow_rect, X(i), Y(i));
end

% Make the image into a texture
stimulus_imageTexture = Screen('MakeTexture', stimulus_window, theImage);
eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);

# grating
grating_gridSize = 256;
orientations = [0:15:180];
pixelsPerPeriod = 33;
plateauCycles = 3;
edgeCycles = 0.25;
contrast = 0.8
windowPointer = stimulus_window;

grating_rect  = [1 1 grating_gridSize grating_gridSize];
[gXoff,gYoff] = pol2cart(deg2rad(135),135);
grating_rect  = CenterRectOnPoint(grating_rect,screenXpixels/2+gXoff,screenYpixels/2+gYoff);

grText     = generate_grating_textures(gridSize,orientations,pixelsPerPeriod,plateauCycles,edgeCycles,stimulus_window,stimulus_screenNumber,contrast);
grText_iTR = generate_grating_textures(gridSize,orientations,pixelsPerPeriod,plateauCycles,edgeCycles,eyeTrack_window,eyeTrack_screenNumber,contrast);

fix_point_rect = [1 1 10 10];
fix_point_rect = CenterRectOnPoint(fix_point_rect, screenXpixels/2, screenYpixels/2);

fix_point_Window = s1*scaler*trackWin_factor;
fix_point_Windowrect = [1 1 fix_point_Window fix_point_Window];
fix_point_Windowrect = CenterRectOnPoint(fix_point_Windowrect,screenXpixels/2, screenYpixels/2);
fix_point_Window = fix_point_Window/2;

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

  gaze_position = gaze_position*nan;
  g = 0;
  
  tr_ind = tr_ind + 1
  
  waiting_for_response = 0;
  tracking_reward = 0;
  on_target = 0;
  text_ind = randi(length(grText_iTR));
  min_target_time = response_wait_min + (response_wait_max - response_wait_min)*rand();
  
  if tr_ind > 1
    Screen('Flip', eyeTrack_window,0);
  end   
  
  % Grey screen
  Screen('FillRect', stimulus_window, grey, rects(:,:,pos));    
  Screen('FillRect', stimulus_window, grey, grating_rect);  
  
  Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));  
  Screen('FillRect', eyeTrack_window, grey, grating_rect);  
  
  greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);  
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
      save(saveSTR,'trial_records');
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
    
    g = g +1;
    gaze_position(:,g) = XY;
    
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) < trackWindow
      on_target = 1;
      reaction_time = toc(wait_fixation_clock);
      eyeTrack_clock = tic();      
      
      Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));                         
      Screen('DrawTexture', eyeTrack_window, grText_iTR(text_ind), [], grating_rect);
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));
      Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
      
      Screen('FillRect', stimulus_window, grey, rects(:,:,pos));                      
      Screen('DrawTexture', stimulus_window, grText(text_ind), [], grating_rect);
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos)); 
      
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
    
    g = g +1;
    gaze_position(:,g) = XY;
    
    if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow      
      on_target = 0;       
      on_target_time   = toc(eyeTrack_clock);
      reward_size_time = toc(eyeTrack_clock);      
      trial_error = 'broke_early';
      if reward_size_time > min_fixation_time
        reward_size_time = reward_scaler*sqrt(reward_size_time);
        Datapixx('SetDoutValues', 1);
        Datapixx('RegWrRd');
        a = tic();      
        while toc(a) < reward_size_time
          # pump juice
        end
        # turn pump off
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        trial_error = 'no_error';        
        break;      
      else
        break;
      end
    end    
   
    # the animals needs to hold fixation for "min_target_time" before responding
    if toc(eyeTrack_clock) > max_fixation_time
      on_target = 0;
      on_target_time   = toc(eyeTrack_clock);
      reward_size_time = toc(eyeTrack_clock);
      reward_size_time = reward_scaler*sqrt(reward_size_time);      
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();      
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');
      trial_error = 'no_error';
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;
      close all;
      sca;
      save(saveSTR,'trial_records');
      break;
    end    
  end  %   
  
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      save(saveSTR,'trial_records');
      break;
  end
 
  % record trial data
  trial_records(tr_ind).stimulus = stimulus_image;
  trial_records(tr_ind).positionX = X(pos);
  trial_records(tr_ind).positionY = Y(pos);
  trial_records(tr_ind).reaction_time  = reaction_time;
  trial_records(tr_ind).on_target_time = on_target_time; 
  trial_records(tr_ind).trial_error = trial_error;
  trial_records(tr_ind).gaze_position = gaze_position;
  
  save(['/home/vpixx/MonkeyRecords/TrialRecords/',animal,'/',saveSTR],'trial_records');
  
  if tr_ind >= max_trs
    is_running = 0;
    sca;
  end

end