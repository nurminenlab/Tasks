% prepare PsychoToolbox
addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
sca;
close all;

# use mouse instead of eye tracker
mouse_track = 1;
debug_on = 0;
save_records = 0;
animal = 'Sansa';
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


  
Trans_mx_shift = [30 -45]; # a manual offset to the translation matrix of the eye tracker calibration. DEF in pixels. 

% task parameters
fix_target_deg       = 1.8;
fix_target_pix       = fix_target_deg*va_in_pixels;
track_win_deg        = 3;
track_win_pix        = track_win_deg*va_in_pixels;
d_target_extra       = 2.5;

wait_fixation        = 1;
rewardConsume_period = 2;
max_fixation_time    = 2;
min_target_time      = 0.2;  
gaze_move_time       = 0.4;
response_wait_time   = gaze_move_time;
max_trs              = 1000;
wrong_target_abort = 1;

gridSize = 100;
contrasts = [0.64];
contrasts_idx = [1]; # workaround for weighted randomization of contrast  
orientations = [0,90];
pix_per_period = 33;
plateau_deg = 1.5;
plateau_pix = plateau_deg*va_in_pixels;
edge_deg = 0.1;
edge_pix = edge_deg*va_in_pixels;  
d_target = (plateau_deg + edge_deg + d_target_extra)*va_in_pixels;
reward_scaler = 0.5;
animal_code = 'MM001';  

thetas_deg = [-45,45,135,225] + 0;
R_deg  = 4;
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
    
    
  % Draw monkey face
  Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
  Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));    
  % Draw fixation window  
  Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
  Datapixx('SetDoutValues', 4);
  stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
  stimulusScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window, greyScreen_eyeTrack_vbl + rewardConsume_period,1);
  Datapixx('RegWrRd');
  

  on_target = 1;
    
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

  Screen('Flip', stimulus_window,0,1);     
  Screen('Flip', eyeTrack_window,0,1);  
  imageArray = Screen('GetImage',stimulus_window);
  is_running = 0;
  sca;
end