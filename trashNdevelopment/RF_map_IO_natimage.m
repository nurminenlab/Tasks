function trial_records = RF_map_IO(debug_on)
  % prepare PsychoToolbox
  addpath('/home/vpixx/Tasks/Functions/');
  AssertOpenGL;
  sca;
  close all;

  # use mouse instead of eye tracker
  mouse_track = 1;

  animal = 'Sansa-';
  saveSTR = [animal,'RF-map-trial_records-',date,'.mat'];
  save_append = 0;
  while exist(saveSTR,'file') == 2
    save_append = save_append + 1;
    saveSTR = [animal,'RF-map-trial_records-',date,'_v',num2str(save_append),'.mat'];
  end
  

  % user defined parameters
  if strcmp(animal,'Tully-')
    scaler               = 0.6;
    small_scaler         = 0.6;
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
    scaler               = 0.4; # KEVIN can you please modify so that we can define the sizes of fixation targets and other stimulus in degrees of visual angle
    small_scaler         = 0.4;
    trackWin_factor      = 2.2;
    %raised cosine variables
    scale_val1 = 200; 
    scale_val2 = 20;    
    %stimulus image scaling
    stimulus_size = tan(1/2) * 47; # KEVIN units in pixels
    stimulus_scaler = 0.8 * stimulus_size;
    stimulus_center = [825 675];
    %stimulus_size = 150;
    image_duration = 0.25;
    waitframes = ceil(image_duration*120);
    %turn on/off fixation mask 
    fill_fixation = 1;
    %turn on/off black white images
    black_white = 0;
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
  
  # 
  expt_info.scaler = scaler;
  expt_info.smallScaler = small_scaler;
  expt_info.trackWinFactor = trackWin_factor;
  expt_info.scaleVal1 = scale_val1;
  expt_info.scaleVal2 = scale_val2;
  expt_info.stimulusScaler = stimulus_scaler;
  expt_info.stimulusCenter = stimulus_center;
  expt_info.stimulusScaler = stimulus_scaler;
  expt_info.waitframes = waitframes;
  expt_info.waitFixation = wait_fixation;
  expt_info.rewardConsumePeriod = rewardConsume_period;
  expt_info.ms = ms;
  expt_info.minTargetTime = min_target_time;
  expt_info.maxTrs = max_trs;
  expt_info.responseWaitMin = response_wait_min;
  expt_info.responseWaitMax = response_wait_max;
  expt_info.gazeMoveTime = gaze_move_time;
  expt_info.maxFixationTime = max_fixation_time;
  expt_info.minFixationTime = min_fixation_time;
  expt_info.rewardScaler = reward_scaler;
  expt_info.FR = FR;
  expt_info.gridSize = gridSize;
  expt_info.fixPointWindowSize = fix_point_Window_size;
  
  
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
  
  %% NEW
  ifi1 = Screen('GetFlipInterval', eyeTrack_window);
  ifi2 = Screen('GetFlipInterval', stimulus_window);
  Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
  
  if(fill_fixation == 1)
    [track_centerX, track_centerY] = WindowCenter(eyeTrack_window);
    fillStimulus = [track_centerX-43, track_centerY-43, track_centerX+43, track_centerY+43];
  else
    fillStimulus = [0 0 0 0];
  end
  
  numImage = 300;
  theImages = cell(1, numImage);
  imageTextures = cell(1, numImage);
  theImageLocations = cell(1, numImage);
  
  stimulus_rect = [stimulus_center(1)-stimulus_size*stimulus_scaler/2, stimulus_center(2)-stimulus_size*stimulus_scaler/2,...
                  stimulus_center(1)+stimulus_size*stimulus_scaler/2, stimulus_center(2)+stimulus_size*stimulus_scaler/2]; 

  
  cd; # remove cds back and forth, sorry maybe my idiosyncratic style
  cd ImagesandRecords/images;
  fls = dir(give path); # KEVIN you can give the full path as an argument
  
  # make textures
  for i = 1:numImage
    % Here we load in an image from file. This one is a image of rabbits that
    % is included with PTB, then shift the images to blur them into the background
    theImageLocations{i} = [add full path, '/', fls(i+2).name; # include fls.path + '/'
    theImages{i} = imread(theImageLocations{i});
    theImages{i} = theImages{i} - 128;

    % Resize and Get the size of the image
    [s1, s2, s3] = size(theImages{i});

    % Here we check if the image is too big to fit on the screen and abort if
    % it is. See ImageRescaleDemo to see how to rescale an image.
    if s1 > screenYpixels || s2 > screenYpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        return;
    end
 
    % Create the cosine window
    [cos_window] = raised_cosine(scale_val1,scale_val2,s2,s1,0,0,'R');        
    
    % Alpha Blend the images with cosine window
    alpha = zeros(s1,s2,3);
    alpha(:,:,1) = cos_window;
    alpha(:,:,2) = cos_window;
    alpha(:,:,3) = cos_window;
    theImages{i} = theImages{i} .* alpha;
    
    theImages{i} = theImages{i} + 128;
    
    if black_white == 1
      theImages{i}(:,:,2) = theImages{i}(:,:,1);
      theImages{i}(:,:,3) = theImages{i}(:,:,1);
    end

    % Make the image into a texture
    imageTextures{i} = Screen('MakeTexture', stimulus_window, theImages{i});
    imageTextures{i} = Screen('MakeTexture', eyeTrack_window, theImages{i});

  end
  
  expt_info.imageLocation = theImageLocations;
    
  counter = randperm(300, 300);  # KEVIN no hard coding
  count = 1;
  trCount = [];
  
  % Load marmoset face fixation target
  stimulus_image = 'face10.jpg';
  theImage = imread(stimulus_image); 
  
  [s1, s2, s3] = size(theImage);
  
  % scale image rectangle
  rect = [0 0 s1*scaler s2*scaler];
  small_rect = [0 0 s1*small_scaler s2*small_scaler];
  
  eyePos_rect = [0 0 5 5];
  trackWindow_rect = [0 0 s1*scaler*trackWin_factor s2*scaler*trackWin_factor];
  trackWindow = s1*scaler*trackWin_factor/2;
  
  idx = 0;
  XBC = [screenXpixels/2-125, screenXpixels/2, screenXpixels/2 + 125];
  YBC = [screenYpixels/2-125, screenYpixels/2, screenYpixels/2 + 125];
  [X,Y] = meshgrid(XBC,YBC);
  
  for i = 1:numel(X);  
    rects(:,:,i) = CenterRectOnPoint(rect, X(i), Y(i));
    small_rects(:,:,i)  = CenterRectOnPoint(small_rect, X(i), Y(i));
    trackWindow_rect(:,:,i) = CenterRectOnPoint(trackWindow_rect, X(i), Y(i));
  end
  
  % Make the image into a texture
  stimulus_imageTexture = Screen('MakeTexture', stimulus_window, theImage);
  eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);
  
  # grating
  # KEVIN please remove code that is not used anymore 
  grating_gridSize = 256;
  orientations = [0:15:180];
  pixelsPerPeriod = 33;
  plateauCycles = 3;
  edgeCycles = 0.25;
  contrast = 0.8;
  windowPointer = stimulus_window;
  
  %grating_rect  = [1 1 grating_gridSize grating_gridSize];
  grating_rect = [1 1 300 300];
  [gXoff,gYoff] = pol2cart(deg2rad(135),135);
  # KEVIN remove unnecessary code  
  grating_rect  = CenterRectOnPoint(grating_rect,screenXpixels/2+gXoff,screenYpixels/2+gYoff);
  
  # KEVIN remove unnecessary code
  grText     = generate_grating_textures(gridSize,orientations,pixelsPerPeriod,plateauCycles,edgeCycles,stimulus_window,stimulus_screenNumber,contrast);
  grText_iTR = generate_grating_textures(gridSize,orientations,pixelsPerPeriod,plateauCycles,edgeCycles,eyeTrack_window,eyeTrack_screenNumber,contrast);
  
  fix_point_rect = [1 1 10 10];
  fix_point_rect = CenterRectOnPoint(fix_point_rect, screenXpixels/2, screenYpixels/2);
  
  fix_point_Window = s1*small_scaler*trackWin_factor;
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
  
  # experiment main loop starts here
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
    
    % Grey screen (stay grey after exiting program)
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
        # Kevin may be style issue, but no cds back and forth please
        cd;
        cd ImagesandRecords/Records;
        save(saveSTR,'expt_info'); # KEVIN edit saveSTR so as to include full path
        cd;
        cd Tasks/FixationCalibrationTraining/PTBTesting;
        break;
      end
      
    endwhile
    
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
      
      
      # this condition checks if the monkey's gaze is on target
      if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) < trackWindow 
        # set digital out to 1 for channel XYZ
        Datapixx('SetDoutValues', 16);
        Datapixx('RegWrRd');
        on_target = 1;
        tic();
        reaction_time = toc(wait_fixation_clock);
        eyeTrack_clock = tic();      
               
        if(count == 300)
          count = 1;
          counter = randperm(300, 300);
        else      
          count = count + 1;
        end
        
        startCounter = count;


        %Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));                         
        Screen('DrawTexture', eyeTrack_window, imageTextures{counter(count)}, [], stimulus_rect, 0);
        Screen('FillOval', eyeTrack_window, [128 128 128], fillStimulus);
        Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], small_rects(:,:,pos));
        Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
       
        
        %Screen('FillRect', stimulus_window, grey, rects(:,:,pos));                      
        Screen('DrawTexture', stimulus_window, imageTextures{counter(count)}, [], stimulus_rect, 0);
        Screen('FillOval', stimulus_window, [128 128 128], fillStimulus);
        Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], small_rects(:,:,pos)); 
        
        
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
      
    end  

    vbl1 = Screen('Flip', eyeTrack_window);
    vbl2 = Screen('Flip', stimulus_window);
    
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
      
      
     
      
      if(count == 300)
        count = 1;
        counter = randperm(300, 300);
      else      
        count = count + 1;
      end
       
      
      %Screen('FillRect', eyeTrack_window, grey, trackWindow_rect(:,:,pos));                         
      Screen('DrawTexture', eyeTrack_window, imageTextures{counter(count)}, [], stimulus_rect, 0);
      %Screen('FillArc', eyeTrack_window, [0 0 0], [], wedge, 5);
      Screen('FillOval', eyeTrack_window, [128 128 128], fillStimulus);
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], small_rects(:,:,pos));      
      Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
      
      %Screen('FillRect', stimulus_window, grey, rects(:,:,pos));                      
      Screen('DrawTexture', stimulus_window, imageTextures{counter(count)}, [], stimulus_rect, 0);
      %Screen('FillArc', stimulus_window, [0 0 0], [], wedge, 5);
      Screen('FillOval', stimulus_window, [128 128 128], fillStimulus);
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], small_rects(:,:,pos));
      
      endCounter = count;
      
      vbl1 = Screen('Flip', eyeTrack_window, vbl1 + (waitframes - 0.5) * ifi1);
      vbl2 = Screen('Flip', stimulus_window, vbl2 + (waitframes - 0.5) * ifi2);
      
      
      if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow      
        # set digital out to 0 for channel XYZ
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        on_target = 0;       
        on_target_time   = toc(eyeTrack_clock);
        
        Screen('FillRect', eyeTrack_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', eyeTrack_window); 
        
        Screen('FillRect', stimulus_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', stimulus_window); 
        
        reward_size_time = toc(eyeTrack_clock);      
        trial_error = 'broke_early';
        if reward_size_time > min_fixation_time
          reward_size_time = reward_scaler*sqrt(reward_size_time);
          # look at this command
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
        # set digital out to 0 for channel XYZ
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        on_target = 0;
        on_target_time   = toc(eyeTrack_clock);
        
        Screen('FillRect', eyeTrack_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', eyeTrack_window); 
        
        Screen('FillRect', stimulus_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', stimulus_window); 
        
        reward_size_time = toc(eyeTrack_clock);# set digital out to 0 for channel XYZ
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
      #Screen('Flip', eyeTrack_window,0,1);
      
      if keyIsDown && KbName(keyCode) == 'q';
        is_running = 0;      
        close all;
        sca;
        cd;
        cd ImagesandRecords/Records;
        save(saveSTR,'expt_info');
        cd;
        cd Tasks/FixationCalibrationTraining/PTBTesting;
        break;
      end    
    end  %   
   
    # KEVIN causes crash if endConter is not initialized
    for i = 1 : (endCounter - startCounter + 1)
      trCount(i) = counter(startCounter + i - 1);
    end
   
    % record trial data
    trial_records(tr_ind).stimulus = stimulus_image;
    trial_records(tr_ind).trCount = trCount;
    %trial_records(tr_ind).counter = counter;
    trial_records(tr_ind).positionX = X(pos);
    trial_records(tr_ind).positionY = Y(pos);
    trial_records(tr_ind).reaction_time  = reaction_time;
    trial_records(tr_ind).on_target_time = on_target_time; 
    trial_records(tr_ind).trial_error = trial_error;
    trial_records(tr_ind).gaze_position = gaze_position; 
  
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca;
      # KEVIN please remove unnecessary cds 
      cd;
      cd ImagesandRecords/Records;
      expt_info.trial_records = trial_records;
      save(saveSTR,'expt_info');
      cd;
      cd Tasks/FixationCalibrationTraining/PTBTesting;
      break;
    end  
    
    
    if tr_ind >= max_trs
      is_running = 0;
      sca;     
    end   

  end
  # last line