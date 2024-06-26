function imageArray = BSD_image_flash(debug_on);

 
  % prepare PsychoToolbox
  addpath('/home/vpixx/Tasks/Functions/');
  
  AssertOpenGL;
  sca;
  close all;

  image_dir = '/home/vpixx/Images/BSD/';
  
  # use mouse instead of eye tracker
  mouse_track = 1;
  save_records = 0;  

  % user defined parameters
    
  distance = 47;
  pix_per_cm = 36.2;
  va_in_pix  = va2pix(distance,pix_per_cm);
  Trans_mx_shift = [30 -45]; # a manual offset to the translation matrix of the eye tracker calibration. DEF in pixels. 
  fill_fixation = 1;    
  black_white = 0;
  
  fixation_target_deg = 1.1;      
  trackWin_deg = 1.6;
 
  stimulus_size_deg = [16];    
  edge_rolloff_deg  = 0.2;
  stimulus_center = [825 675]; # change this so as to be defined in polar coordinates MAYBE LATER
  
  image_duration = 0.3;
  blank_duration = 0.5;
  blank = 0;    
  waitframes = ceil(image_duration*120);    
  waitframes2 = ceil(blank_duration*120);        
  
  wait_fixation        = 0.75;
  rewardConsume_period = 2;
  ms                   = 10;    
  max_trs              = 1;    
  max_fixation_duration    = 24;
  min_fixation_duration    = 0.1;
  reward_scaler = 0.7;
  FR = 120;    
  fix_point_Window_size = 100;
  trackMarkerColor = [255,0,0];
  gaze_position = nan*ones(2,FR*ceil((wait_fixation+max_fixation_duration)));
  TTLwidth = 0.15;
   
  # raised cosine mask parameters  
  mask_grid = 962;
  plateau_pix = stimulus_size_deg*va_in_pix;
  edge_pix = edge_rolloff_deg*va_in_pix;
  
  # expt params
  expt_info.distance   = distance;
  expt_info.pix_per_cm = pix_per_cm;
  expt_info.va_in_pix = va_in_pix;
  expt_info.Trans_mx_shift = Trans_mx_shift;  
  expt_info.fill_fixation = 1;
  expt_info.black_white = 0;
  
  expt_info.fixation_target_deg = fixation_target_deg;
  expt_info.trackWin_deg = trackWin_deg;
  
  expt_info.stimulus_size_deg = stimulus_size_deg;
  expt_info.edge_rolloff_deg = edge_rolloff_deg;  
  expt_info.stimulus_center = stimulus_center;
  
  expt_info.image_duration = image_duration;
  expt_info.blank_duration = blank_duration;  
  expt_info.waitframes = waitframes;
  expt_info.wait_fixation = wait_fixation;
  expt_info.rewardConsume_period  = rewardConsume_period;  
  expt_info.max_fixation_duration = max_fixation_duration;
  expt_info.min_fixation_duration = min_fixation_duration;
  
  expt_info.max_trs = max_trs;  
  expt_info.reward_scaler = reward_scaler;
  expt_info.FR = FR;
  expt_info.mask_grid = mask_grid;
  
  if mouse_track
    XY = ones(2,1)*nan;
  end  
  
  Datapixx('Open');
  adcRate = 1e3;
  baseBuffAddr = 0;
  minStreamFrames = 15;
   
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

  ifi1 = Screen('GetFlipInterval', eyeTrack_window);
  ifi2 = Screen('GetFlipInterval', stimulus_window);
  Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
  Screen('BlendFunction', stimulus_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
  if(fill_fixation == 1)
    [track_centerX, track_centerY] = WindowCenter(eyeTrack_window);
    fillStimulus = [track_centerX-(trackWin_deg*va_in_pix)/2, ...
                      track_centerY-(trackWin_deg*va_in_pix)/2, ...
                      track_centerX+(trackWin_deg*va_in_pix)/2, ...
                      track_centerY+(trackWin_deg*va_in_pix)/2];
  else
    fillStimulus = [0 0 0 0];
  end
  
  alpha_mask = ones(mask_grid,mask_grid,2)*grey;
  mask_text_eyeTrack = NaN*ones(1,length(stimulus_size_deg));  
  mask_text_stimulus = NaN*ones(1,length(stimulus_size_deg));
  for i = 1:length(stimulus_size_deg)
    [cos_window] = 255 - 255*raised_cosine(plateau_pix(i),edge_pix,mask_grid,mask_grid,0,0,'R');  
    alpha_mask(:,:,2) = cos_window;    
    mask_text_eyeTrack(i) = Screen('MakeTexture', eyeTrack_window, alpha_mask);
    mask_text_stimulus(i) = Screen('MakeTexture', stimulus_window, alpha_mask);
  end
  
  mask_rect = [0 0 mask_grid mask_grid];
  mask_rect = CenterRectOnPoint(mask_rect, stimulus_center(1), stimulus_center(2));
  
  specificImage = [40,13];
  #specificImage = [1 3 5 9 11 13];
  numImage = length(specificImage);  
  fls = dir(image_dir);
  
  all_rects = ones(numImage,4)*NaN;
  rect_type = ones(numImage,1)*NaN;
  # make textures
  for i = 1:numImage
    % Here we load in an image from file. This one is a image of rabbits that
    % is included with PTB, then shift the images to blur them into the background
    theImageLocations{i} = [fls(specificImage(i)+2).folder,'/',fls(specificImage(i)+2).name];    
    theImages{i} = imread(theImageLocations{i});
    theImages{i} = theImages{i};

    % Resize and Get the size of the image
    [s1, s2, s3] = size(theImages{i});
    
    if s1 > s2
      rect_type(i) = 1;
    elseif s1 < s2
      rect_type(i) = 2;          
    endif
    
    all_rects(i,:) = [0 0 2*s1 2*s2];
    
    if black_white == 1      
      theImages{i} = mean(theImages{i},3);
    end

    % Make the image into a texture
    imageTextures_stimulus{i} = Screen('MakeTexture', stimulus_window, theImages{i});
    imageTextures_eyeTrack{i} = Screen('MakeTexture', eyeTrack_window, theImages{i});

  end  

  #all_rects = unique(all_rects,'rows');  
  all_rects = round(all_rects);
  
  # translate to correct location
  for i = 1:size(all_rects,1)
    all_rects(i,:) = CenterRectOnPoint(all_rects(i,:),stimulus_center(1), stimulus_center(2));
  endfor 
  
  expt_info.imageLocation = theImageLocations;
  
  Screen('TextSize', eyeTrack_window, 50);
  Screen('TextFont', eyeTrack_window, 'Courier');
  DrawFormattedText(eyeTrack_window, "Ready", screenXpixels/2 - 70, screenYpixels/2, [256 0 0]);
  Screen('Flip', eyeTrack_window);    
  
  counter = randperm(numImage, numImage); 
  count = 0;
  trCount = [];
  
  
  imgScaleOrder = randperm(length(stimulus_size_deg),length(stimulus_size_deg));
  scaleCount = 1;     
  
  % Load marmoset face fixation target
  stimulus_image = 'face8.jpg';
  theImage = imread(stimulus_image);
  [s1, s2, s3] = size(theImage);
  
  % scale image rectangle
  rect = [0 0 va_in_pix*fixation_target_deg va_in_pix*fixation_target_deg];  
  
  eyePos_rect = [0 0 5 5];
  trackWindow_rect = [0 0 trackWin_deg*va_in_pix trackWin_deg*va_in_pix]; 
  trackWindow = va_in_pix*trackWin_deg/2;
  
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
  
  fix_point_Window = va_in_pix*trackWin_deg;
  fix_point_Windowrect = [1 1 fix_point_Window fix_point_Window];
  fix_point_Windowrect = CenterRectOnPoint(fix_point_Windowrect,screenXpixels/2, screenYpixels/2);
  fix_point_Window = fix_point_Window/2;
  
  tr_ind = 0;
  trial_records = struct();
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
  
  try
    # experiment main loop starts here
    while is_running
    
      if(scaleCount == length(stimulus_size_deg))
          scaleCount = 1;
          imgScaleOrder = randperm(length(stimulus_size_deg), length(stimulus_size_deg));
      else      
          scaleCount = scaleCount + 1;
      end
    
      gaze_position = gaze_position*nan;
      g = 0;
    
      tr_ind = tr_ind + 1
    
      waiting_for_response = 0;
      tracking_reward = 0;
      on_target = 0;      
      even = 0;
      
      if tr_ind > 1
        Screen('Flip', eyeTrack_window,0);
      end   
    
      % Grey screen (stay grey after exiting program)
      Screen('FillRect', stimulus_window, grey, rects(:,:,pos));    
      Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));  
    
      greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);  
      greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window);    
      
    
      % Draw monkey face
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));     
      
      % Draw fixation window
      Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);      
      stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period);
      stimulusScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window, greyScreen_eyeTrack_vbl + rewardConsume_period,1);      
      
          
      # set digital out to 1 for channel XYZ when the monkey acquires fixation                   
      trImage = [];
      trImage_idx = 1;
      trScale = [];
      trScale_idx = 1;
      
      if(count == numImage)
        count = 1;
        counter = randperm(numImage, numImage);
      else      
        count = count + 1;
      end      

      stimulus_rect = all_rects(rect_type(counter(count)),:);          
      startCounter = count;

      # first stimulus after the monkey acquired fixation
      Screen('DrawTexture', eyeTrack_window, imageTextures_eyeTrack{counter(count)}, [], stimulus_rect, 0);
      Screen('FillOval', eyeTrack_window, [128 128 128], fillStimulus);
      Screen('DrawTexture', eyeTrack_window, mask_text_eyeTrack(imgScaleOrder(scaleCount)), [], mask_rect, 0);
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));
      Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
   
      Screen('DrawTexture', stimulus_window, imageTextures_stimulus{counter(count)}, [], stimulus_rect, 0);
      Screen('FillOval', stimulus_window, [128 128 128], fillStimulus);
      Screen('DrawTexture', stimulus_window, mask_text_stimulus(imgScaleOrder(scaleCount)), [], mask_rect, 0);
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos)); 
    
      vbl1 = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);     
      vbl2 = Screen('Flip', eyeTrack_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
        
      imageArray = Screen('GetImage',stimulus_window);     
      sca;      
      break;
    end 
  
catch 
  lasterr
  sca;  
end
end

  # last line