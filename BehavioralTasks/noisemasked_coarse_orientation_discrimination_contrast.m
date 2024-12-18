% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;
contrast = 0;

save_default_options('-mat7-binary');

# use mouse instead of eye tracker

mouse_track = 0;
debug_on = 0;
save_records = 1;

animal = 'Wolfjaw';

records_folder = '/home/vpixx/MonkeyRecords/TrialRecords/';
saveSTR = [records_folder,animal,'/','noisemasked_coarse_orientation_discrimination-trial_records-',date,'.mat'];
save_append = 0;
while exist(saveSTR,'file') == 2
  save_append = save_append + 1;
  saveSTR = [records_folder,animal,'/','noisemasked_coarse_orientation_discrimination-trial_records-',date,'_v',num2str(save_append),'.mat'];
end

distance = 47;
pix_per_cm = 36.2; 
va_in_pixels = va2pix(distance,pix_per_cm);

if strcmp(animal,'Sansa')
  
  Trans_mx_shift = [30 -15];# a manual offset to the translation matrix of the eye tracker calibration. DEF in pixels. 
  
  % task parameters
  fix_target_deg       = 0.74;
  fix_target_pix       = fix_target_deg*va_in_pixels;
  track_win_deg        = 1.5;
  track_win_pix        = track_win_deg*va_in_pixels;
  d_target_extra       = 2.5;
  R_deg  = 4;
  rotate_deg           = 0; # goal to use something close to 30, I hope
  
  wait_fixation        = 1;
  rewardConsume_period = 2;
  max_fixation_time    = 2;
  min_target_time      = 0.2;  
  gaze_move_time       = 0.4;
  response_wait_time   = gaze_move_time;
  max_reward_time      = 0.75;
  max_trs              = 10000;
  wrong_target_abort = 1;
  
  gridSize = 128;
  #contrasts = [0.01 0.02 0.04 0.08 0.16 0.32];
  #contrasts = [0.01 0.02 0.04 0.08 0.16 0.32];
  contrasts = [0.04 0.08 0.16 0.32];
  contrasts_idx = [1 2 3 4]; # workaround for weighted randomization of contrast  
  #contrasts_idx = [1 1 2 3 4 5 6];
  orientations = [0,90];
  pix_per_period = 33;
  plateau_deg = 3;
  plateau_pix = plateau_deg*va_in_pixels;
  edge_deg = 0.1;
  edge_pix = edge_deg*va_in_pixels;  
  d_target = (plateau_deg + edge_deg + d_target_extra)*va_in_pixels;
  reward_scaler = 0.53;
  animal_code = 'MM001';  
  wrong_target_timeout = 0.5;
  no_fixation_timeout = 1;
  
  driveAP = nan;
  driveML = nan;
  sex = 'F';
  DOB = nan;
  handler = 'SC';
  experimenters = 'SC'; 
  weight = 463;   
  penetration_wait_time = nan;
  penetration_time = nan;
  
  
elseif strcmp(animal,'Wolfjaw')
  
  Trans_mx_shift       = [0 0]; # a manual offset to the translation matrix of the eye tracker calibration. DEF in pixels. 
  
  % task parameters
  fix_target_deg       = 1; # fixation target size in degrees
  fix_target_pix       = fix_target_deg*va_in_pixels;
  track_win_deg        = 2.25; # blue ring around central fixation target
  track_win_pix        = track_win_deg*va_in_pixels;
  d_target_extra       = 2.5; # acceptance window for other targets
  R_deg  = 4;
  rotate_deg           = 0;
  
  wait_fixation        = 1; # wait fixation time until next trial offered
  rewardConsume_period = 1; # time between trials
  max_fixation_time    = 1.5; # length of time that the monkey can look at the larger target
  min_target_time      = 0.13; # min fixation time to central target before giving answer; increase to 0.2 gradually
  gaze_move_time       = 0.4; 
  response_wait_time   = gaze_move_time;
  max_reward_time      = 1.45;
  max_trs              = 10000;
  wrong_target_abort = 1;
  
  gridSize = 128;
  #contrasts = [0.02 0.04 0.08 0.16 0.32];
  contrasts = [0.04 0.08 0.16 0.32];
  #contrasts_idx = [1,1,2,2,3,3,4,5];# workaround for weighted randomization of contrast
  contrasts_idx = [1 1 2 3 4];
  orientations = [0,90];
  pix_per_period = 33;
  plateau_deg = 3;
  plateau_pix = plateau_deg*va_in_pixels;
  edge_deg = 0.1;
  edge_pix = edge_deg*va_in_pixels;  
  d_target = (plateau_deg + edge_deg + d_target_extra)*va_in_pixels;
  reward_scaler = 0.6;
  animal_code = 'MM004';
  wrong_target_timeout = 0.5;
  no_fixation_timeout = 1;
  
  driveAP = nan;
  driveML = nan;
  sex = 'M';
  DOB = nan;
  handler = 'SC';
  experimenters = 'SCMM'; 
  weight = 429;   
  penetration_wait_time = nan;
  penetration_time = nan;
  
else

  display('ERROR: No such animal')  
  return;

endif

thetas_deg = [-45,45,135,225] + rotate_deg;
thetas = deg2rad(thetas_deg);
R      = va_in_pixels*R_deg;
ms     = 10;

expt_info.project              = 'inhibitory-neuron-opto';
expt_info.experiment           = 'masked_discrimination';
expt_info.date                 = date();
expt_info.animal               = animal;
expt_info.animal_code          = animal_code;
expt_info.fix_target_deg       = fix_target_deg;
expt_info.track_win_deg        = track_win_deg;
expt_info.d_target_extra       = d_target_extra;
expt_info.wait_fixation        = wait_fixation;
expt_info.rewardConsume_period = rewardConsume_period;
expt_info.max_fixation_time    = max_fixation_time;
expt_info.min_target_time      = min_target_time;
expt_info.gaze_move_time       = gaze_move_time;
expt_info.response_wait_time   = response_wait_time;
expt_info.wrong_target_abort   = wrong_target_abort;
expt_info.contrasts            = contrasts;
expt_info.contrasts_idx        = contrasts_idx;
expt_info.orientations         = orientations;
expt_info.pix_per_period       = pix_per_period
expt_info.plateau_deg          = plateau_deg;
expt_info.edge_deg             = edge_deg;
expt_info.reward_scaler        = reward_scaler;
expt_info.wrong_target_timeout = wrong_target_timeout;
expt_info.no_fixation_timeout  = no_fixation_timeout;
expt_info.max_reward_time      = max_reward_time;
expt_info.rotate_deg           = rotate_deg;

fix_point_Window_size = track_win_pix;
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
  Trans_mx = [bx(1)+Trans_mx_shift(1), by(1)+Trans_mx_shift(2)]';
else
  Scale_mx = eye(2);
  Trans_mx = [0 0];
end

expt_info.Scale_mx = Scale_mx;
expt_info.Trans_mx = Trans_mx;

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
Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('BlendFunction', stimulus_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the size of the on screen window, these are the same for both screens
[screenXpixels, screenYpixels] = Screen('WindowSize',stimulus_window);

% Load marmoset face
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);

% scale image rectangle
rect = [0 0 fix_target_pix fix_target_pix];
eyePos_rect = [0 0 5 5];
trackWindow_rect = [0 0 track_win_pix track_win_pix];
trackWindow = track_win_pix;

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

target_Windowrect = [1 1 gridSize,gridSize];
target_Windowrect = CenterRectOnPoint(target_Windowrect,screenXpixels/2, screenYpixels/2);

target_MarkRect = [0, 0, d_target, d_target];

txt_rect = [1 1 gridSize gridSize];
txt_rect = CenterRectOnPoint(txt_rect, screenXpixels/2, screenYpixels/2);
txt_rects  = nan * ones(4,length(thetas)*length(R));

[Xoff, Yoff] = pol2cart(thetas,R);
for i = 1:length(Xoff)
  txt_rects(:,i) = round(CenterRectOnPoint(txt_rect, screenXpixels/2+Xoff(i), screenYpixels/2+Yoff(i)));
end

# generate gratings
grText_all     = generate_grating_textures(gridSize,orientations,pix_per_period,stimulus_window,stimulus_screenNumber);
grText_all_iTrack = generate_grating_textures(gridSize,orientations,pix_per_period,eyeTrack_window,eyeTrack_screenNumber);

# generate cosine masks 
alpha_mask = ones(gridSize,gridSize,2)*grey;
mask_text_eyeTrack = NaN*ones(1,length(contrasts));  
mask_text_stimulus = NaN*ones(1,length(contrasts));
for i = 1:length(contrasts)
    [cos_window] = 255 - contrasts(i)*255*raised_cosine(plateau_pix,edge_pix,gridSize,gridSize,0,0,'R');  
    alpha_mask(:,:,2) = cos_window;    
    mask_text_eyeTrack(i) = Screen('MakeTexture', eyeTrack_window, alpha_mask);
    mask_text_stimulus(i) = Screen('MakeTexture', stimulus_window, alpha_mask);
end

tr_ind = 0;
trial_records = struct();
is_running = 1;
conditions = ['5'];
pos = 5;
stim = 'center';
trial_contrast = 1;

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
  #min_target_time = response_wait_min + (response_wait_max - response_wait_min)*rand();
  
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
  
  target_MarkRect = round(CenterRectOnPoint(target_MarkRect, screenXpixels/2+Xoff(target_pos), screenYpixels/2+Yoff(target_pos)));
  
  not_target_Xoff = Xoff;
  not_target_Xoff(target_pos) = [];  
  not_target_pos_X = screenXpixels/2+not_target_Xoff;
  
  not_target_Yoff = Yoff;
  not_target_Yoff(target_pos) = [];
  not_target_pos_Y = screenYpixels/2+not_target_Yoff;
  
  grText = repmat(grText_all(2),1,4);
  grText(:,target_pos) = grText_all(1);
  grText_iTrack = repmat(grText_all_iTrack(2),1,4);
  grText_iTrack(:,target_pos) = grText_all_iTrack(1);
  
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
      if save_records
        expt_info.trial_records = trial_records;
        save(saveSTR,'expt_info');
      endif      
      break;
    end
  end  %
    
  % Draw monkey face
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));    
  % Draw fixation window  
  Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
  Datapixx('SetDoutValues', 4);
  stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
  stimulusScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window, greyScreen_eyeTrack_vbl + rewardConsume_period,1);
  Datapixx('RegWrRd');
  
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
      fixation_reaction_time = toc(wait_fixation_clock);
      eyeTrack_clock = tic();      
      
      contrast = contrasts(contrasts_idx(trial_contrast));
      #Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));           
      #Screen('FillOval', eyeTrack_window,[0 0 255],fix_point_rect);       
      Screen('DrawTextures', eyeTrack_window, grText_iTrack,[],txt_rects);
      Screen('DrawTextures', eyeTrack_window, repmat(mask_text_eyeTrack(contrasts_idx(trial_contrast)),1,4),[],txt_rects);
      Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3);      
      
      #Screen('FillRect', stimulus_window, grey, rects(:,:,pos));           
      #Screen('FillOval', stimulus_window,[0 0 255],fix_point_rect);
      Screen('DrawTextures', stimulus_window, grText,[],txt_rects);
      Screen('DrawTextures', stimulus_window, repmat(mask_text_stimulus(contrasts_idx(trial_contrast)),1,4),[],txt_rects);
      
      Datapixx('SetDoutValues', 20);          
      Screen('Flip', stimulus_window,0,1);     
      Screen('Flip', eyeTrack_window,0,1);
      Datapixx('RegWrRd');
      break;
    end
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    reaction_time = nan;
    fixation_reaction_time = nan;
    trial_error = 'no_fixation';
    rewardConsume_period = no_fixation_timeout;
    on_target_time = nan;
    selected_pos_X = nan;
    selected_pos_Y = nan;
    reaction_time = nan;
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
      trial_error = 'broke_fixation';
      rewardConsume_period = no_fixation_timeout;
      on_target_time = toc(eyeTrack_clock);
      selected_pos_X = nan;
      selected_pos_Y = nan;
      reaction_time = nan;
      break;      
    end           
   
    # the animals needs to hold fixation for "min_target_time" before responding
    if toc(eyeTrack_clock) > min_target_time
      on_target = 0;
      waiting_for_response = 1;
      on_target_time = toc(eyeTrack_clock);
      trial_error = 'no_error';
      Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos)); # blank fixation target?
      #Screen('DrawTextures', eyeTrack_window, grText_iTrack,[],txt_rects);
      Screen('FrameOval',eyeTrack_window,[0 0 255],target_MarkRect,3,3);
      
      Screen('FillRect', stimulus_window, grey, rects(:,:,pos)); # blank fixation target?
      #Screen('DrawTextures', stimulus_window, grText,[],txt_rects);
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      response_time_clock = tic();
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      if save_records
        expt_info.trial_records = trial_records;
        save(saveSTR,'expt_info');
      endif
      
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
      trial_error = 'no_response';
      rewardConsume_period = no_fixation_timeout;
      selected_pos_X = nan;
      selected_pos_Y = nan;
      reaction_time = nan;
      break;
    end
    
    # target hit
    if sqrt((XY(1) - target_pos_X)^2 + (XY(2) - target_pos_Y)^2) < d_target/2
      Screen('FillRect', eyeTrack_window, grey, txt_rects);
      Screen('FillRect', stimulus_window, grey, txt_rects);
      % Draw monkey face
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], txt_rects(:,target_pos));  
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], txt_rects(:,target_pos));
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      reward_time_clock = tic();
      tracking_reward = 1;
      
      trial_error = 'hit'; #no need to setup RewardConsume_period
      selected_pos_X = target_pos_X;
      selected_pos_Y = target_pos_Y;
      reaction_time = toc(response_time_clock);
      
      trial_contrast = trial_contrast + 1;
      if trial_contrast > length(contrasts_idx)
        trial_contrast = 1;
        contrasts_idx = contrasts_idx(randperm(length(contrasts_idx)));        
      endif
      break;
    end
    
    if wrong_target_abort
      if sqrt((XY(1) - not_target_pos_X(1))^2 + (XY(2) - not_target_pos_Y(1))^2) < d_target/2
        waiting_for_response = 0;
        
        trial_error = 'wrong_target';
        rewardConsume_period = wrong_target_timeout;
        selected_pos_X = not_target_pos_X(1);
        selected_pos_Y = not_target_pos_Y(1);
        reaction_time = toc(response_time_clock);
        
        trial_contrast = trial_contrast + 1;
        if trial_contrast > length(contrasts_idx)
          trial_contrast = 1;
          contrasts_idx = contrasts_idx(randperm(length(contrasts_idx)));        
        endif
        break;
      end
      
      if sqrt((XY(1) - not_target_pos_X(2))^2 + (XY(2) - not_target_pos_Y(2))^2) < d_target/2
        waiting_for_response = 0;
        
        trial_error = 'wrong_target';
        rewardConsume_period = wrong_target_timeout;
        selected_pos_X = not_target_pos_X(2);
        selected_pos_Y = not_target_pos_Y(2);
        reaction_time = toc(response_time_clock);
        
        trial_contrast = trial_contrast + 1;
        if trial_contrast > length(contrasts_idx)
          trial_contrast = 1;
          contrasts_idx = contrasts_idx(randperm(length(contrasts_idx)));        
        endif
        break;
      end
      
      if sqrt((XY(1) - not_target_pos_X(3))^2 + (XY(2) - not_target_pos_Y(3))^2) < d_target/2
        waiting_for_response = 0;
        
        trial_error = 'wrong_target';
        rewardConsume_period = wrong_target_timeout;
        selected_pos_X = not_target_pos_X(3);
        selected_pos_Y = not_target_pos_Y(3);
        reaction_time = toc(response_time_clock);
        
        trial_contrast = trial_contrast + 1;
        if trial_contrast > length(contrasts_idx)
          trial_contrast = 1;
          contrasts_idx = contrasts_idx(randperm(length(contrasts_idx)));        
        endif
        break;
      end
    end
  
    # max length of trial
    if toc(response_time_clock) > response_wait_time      
      waiting_for_response = 0; 
      trial_error = 'no_response';
      rewardConsume_period = no_fixation_timeout;
      selected_pos_X = not_target_pos_X(2);
      selected_pos_Y = not_target_pos_Y(2);     
      reaction_time = nan;
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      if save_records
        expt_info.trial_records = trial_records;
        save(saveSTR,'expt_info');
      endif      
      break;
    end    
  end  %

  while tracking_reward;
    
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
   
    if sqrt((XY(1) - target_pos_X)^2 + (XY(2) - target_pos_Y)^2) > 200
      Screen('FillRect', eyeTrack_window, grey);
      Screen('FillRect', stimulus_window, grey);  
      Screen('Flip', stimulus_window, 0);     
      Screen('Flip', eyeTrack_window, 0,1);
      tracking_reward = 0;
      
      reward_size_time = reward_scaler*sqrt((toc(reward_time_clock)));
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();      
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');
      rewardConsume_period = no_fixation_timeout; #To give sometime to consume      
      break;
    end
    
    # max length of trial
    if toc(reward_time_clock) > max_reward_time
      reward_size_time = 0.4*sqrt((toc(reward_time_clock)));
      Datapixx('SetDoutValues', 1);
      Datapixx('RegWrRd');
      a = tic();      
      while toc(a) < reward_size_time
        # pump juice
      end
      # turn pump off
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');
      rewardConsume_period = no_fixation_timeout; #To give sometime to consume   
      tracking_reward = 0;      
      break;
    end 
      
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
    Screen('Flip', eyeTrack_window,0,1);
    
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      if save_records
        expt_info.trial_records = trial_records;
        save(saveSTR,'expt_info');
      endif      
      break;
    end    
  end  %  
  
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      if save_records
        expt_info.trial_records = trial_records;
        save(saveSTR,'expt_info');
      endif      
      break;
  end
 
  % record trial data
  trial_records(tr_ind).stimulus = stimulus_image;
  trial_records(tr_ind).positionX = X(pos);
  trial_records(tr_ind).positionY = Y(pos);
  trial_records(tr_ind).fixation_reaction_time  = fixation_reaction_time;
  trial_records(tr_ind).reaction_time  = reaction_time;
  trial_records(tr_ind).on_target_time = on_target_time; 
  trial_records(tr_ind).trial_error  = trial_error;
  trial_records(tr_ind).target_pos_X = target_pos_X;
  trial_records(tr_ind).target_pos_Y = target_pos_Y;
  trial_records(tr_ind).selected_pos_X = selected_pos_X;
  trial_records(tr_ind).selected_pos_Y = selected_pos_Y;
  trial_records(tr_ind).contrast = contrast;
  trial_records(tr_ind).rewardConsume_period = rewardConsume_period;
  

end