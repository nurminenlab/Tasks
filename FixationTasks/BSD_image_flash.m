function BSD_image_flash(debug_on);
  % prepare PsychoToolbox
  addpath('/home/vpixx/Tasks/Functions/');
  
  AssertOpenGL;
  sca;
  close all;

  image_dir = '/home/vpixx/Images/BSD/';
  
  # use mouse instead of eye tracker
  mouse_track = 0;
  save_records = 0;
  
  animal = 'Sansa';
  saveSTR = ['/home/vpixx/MonkeyRecords/TrialRecords/', animal,'/','BSD-image-flash-trial_records-',date,'.mat'];
  save_append = 0;
  while exist(saveSTR,'file') == 2;
    save_append = save_append + 1;
    saveSTR = ['/home/vpixx/MonkeyRecords/TrialRecords/', animal,'/','BSD-image-flash-trial_records-',date,'_v',num2str(save_append),'.mat'];
  end  

  % user defined parameters
  if strcmp(animal,'Wolfjaw-')    
    # under development
  elseif strcmp(animal,'Sansa')
    distance = 47;
    pix_per_cm = 36.2;
    va_in_pix  = va2pix(distance,pix_per_cm);
    fixation_target_deg = 2;      
    trackWin_deg = 3;
   
    stimulus_size_deg = [4,16];    
    edge_rolloff_deg  = 0.2;
    stimulus_center = [825 675]; # change this so as to be defined in polar coordinates MAYBE LATER
    
    image_duration = 0.3;
    blank_duration = 1;
    blank = 0;    
    waitframes = ceil(image_duration*120);    
    waitframes2 = ceil(blank_duration*120);    
    fill_fixation = 1;    
    black_white = 0;
    
    wait_fixation        = 0.75;
    rewardConsume_period = 2;
    ms                   = 10;    
    max_trs              = 10000;    
    max_fixation_duration    = 24;
    min_fixation_duration    = 0.1;
    reward_scaler = 0.9;
    FR = 120;    
    fix_point_Window_size = 100;
    trackMarkerColor = [255,0,0];
    gaze_position = nan*ones(2,FR*ceil((wait_fixation+max_fixation_duration)));
  else
    fprintf('Strange animal \n');
  end   
  
  # raised cosine mask parameters  
  mask_grid = 962;
  plateau_pix = stimulus_size_deg*va_in_pix;
  edge_pix = edge_rolloff_deg*va_in_pix;
  
  expt_info.fixation_target_deg = fixation_target_deg;
  expt_info.stimulus_size_deg = stimulus_size_deg;
  expt_info.edge_rolloff_deg = edge_rolloff_deg;
  expt_info.va_in_pix = va_in_pix;
  expt_info.trackWin_deg = trackWin_deg;  
  expt_info.stimulus_center = stimulus_center;
  expt_info.waitframes = waitframes;
  expt_info.wait_fixation = wait_fixation;
  expt_info.rewardConsume_period = rewardConsume_period;  
  expt_info.max_trs = max_trs;  
  expt_info.max_fixation_duration = max_fixation_duration;
  expt_info.min_fixation_duration = min_fixation_duration;
  expt_info.reward_scaler = reward_scaler;
  expt_info.FR = FR;
  
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
  
  specificImage = [1 3 5 9 11 13];
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

  all_rects = unique(all_rects,'rows');  
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
  count = 1;
  trCount = [];
  endCounter = 0;
  
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
          if save_records;
            save(saveSTR,'expt_info'); 
          endif
          
          break;
        end
      
      endwhile
    
      % Draw monkey face
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));  
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));     
      
      % Draw fixation window
      Screen('FrameOval',eyeTrack_window, [0 0 255], trackWindow_rect(:,:,pos), 3,3);
      Datapixx('SetDoutValues', 4);
      stimulusScreen_stimulus_vbl = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period);
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
      
        g = g +1;
        gaze_position(:,g) = XY;      
      
        # this condition checks if the monkey's gaze is on target
        if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) < trackWindow 
          
          # set digital out to 1 for channel XYZ when the monkey acquires fixation
          Datapixx('SetDoutValues', 20);
          Datapixx('RegWrRd');
          on_target = 1;
          tic();
          reaction_time = toc(wait_fixation_clock);
          eyeTrack_clock = tic();      
               
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
          Screen('DrawTexture', eyeTrack_window, mask_text_eyeTrack(scaleCount), [], mask_rect, 0);
          Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));
          Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
       
          Screen('DrawTexture', stimulus_window, imageTextures_stimulus{counter(count)}, [], stimulus_rect, 0);
          Screen('FillOval', stimulus_window, [128 128 128], fillStimulus);
          Screen('DrawTexture', stimulus_window, mask_text_stimulus(scaleCount), [], mask_rect, 0);
          Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos)); 
        
          vbl1 = Screen('Flip', stimulus_window, greyScreen_stimulus_vbl + rewardConsume_period,1);     
          vbl2 = Screen('Flip', eyeTrack_window, greyScreen_stimulus_vbl + rewardConsume_period,1);
          # send TTL
          now = tic();
          Datapixx('SetDoutValues', 84);
          Datapixx('RegWrRd');
          while toc(now) < 0.02
            # do nothing
          end
          Datapixx('SetDoutValues', 20);          
          Datapixx('RegWrRd');
          break;
        end
      
        eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
        Screen('FillOval', eyeTrack_window, trackMarkerColor, eyePos_rect);
        Screen('Flip', eyeTrack_window,0,1);
      
        reaction_time = nan;
        trial_error = 'no_fixation';
        on_target_time = nan;
      
      end   
      
    
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
      
        if(count == numImage)
          count = 1;
          counter = randperm(numImage, numImage);
        else      
          count = count + 1;
        end            
        
        even = even + 1;
        if (blank == 1) && (mod(even, 2) == 1)
          stimulus_rect = [0 0 0 0];
          display(even)
          Screen('FillRect', stimulus_window, grey);    
          Screen('FillRect', eyeTrack_window, grey);         
      
          endCounter = count;
      
          vbl1 = Screen('Flip', eyeTrack_window, vbl1 + (waitframes - 0.5) * ifi1);
          vbl2 = Screen('Flip', stimulus_window, vbl2 + (waitframes - 0.5) * ifi2);
        else
          stimulus_rect = all_rects(rect_type(counter(count)),:);         
      
          Screen('DrawTexture', eyeTrack_window, imageTextures_eyeTrack{counter(count)}, [], stimulus_rect, 0);
          Screen('FillOval', eyeTrack_window, [128 128 128], fillStimulus);
          Screen('DrawTexture', eyeTrack_window, mask_text_eyeTrack(scaleCount), [], mask_rect, 0);
          Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));      
          Screen('FrameOval',eyeTrack_window,[0 0 255],fix_point_Windowrect,3,3); # to mark fixation window
      
          Screen('DrawTexture', stimulus_window, imageTextures_stimulus{counter(count)}, [], stimulus_rect, 0);
          Screen('FillOval', stimulus_window, [128 128 128], fillStimulus);
          Screen('DrawTexture', stimulus_window, mask_text_stimulus(scaleCount), [], mask_rect, 0);
          Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));
      
          endCounter = count;
          
          vbl1 = Screen('Flip', eyeTrack_window, vbl1 + (waitframes - 0.5) * ifi1);
          vbl2 = Screen('Flip', stimulus_window, vbl2 + (waitframes - 0.5) * ifi2);
          # send TTL
          now = tic();
          Datapixx('SetDoutValues', 84);
          Datapixx('RegWrRd');
          while toc(now) < 0.02
            # do nothing
          end
          Datapixx('SetDoutValues', 20);
          Datapixx('RegWrRd');
        endif    
        
        if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow      
          # set digital out to 0 for channel XYZ
          Datapixx('SetDoutValues', 0);
          Datapixx('RegWrRd');
          on_target = 0;       
          on_target_time   = toc(eyeTrack_clock);
        
          Screen('FillRect', eyeTrack_window, grey, [0 0 screenXpixels screenYpixels], 0);  
          Screen('FillRect', stimulus_window, grey, [0 0 screenXpixels screenYpixels], 0);
          
          Screen('Flip', eyeTrack_window);              
          Screen('Flip', stimulus_window); 
        
          reward_size_time = toc(eyeTrack_clock);      
          trial_error = 'broke_early';
          if reward_size_time > min_fixation_duration
            reward_size_time = reward_scaler*sqrt(reward_size_time);
            # turn
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
        if toc(eyeTrack_clock) > max_fixation_duration
          # set digital out to 0 for channel XYZ
          Datapixx('SetDoutValues', 0);
          Datapixx('RegWrRd');
          on_target = 0;
          on_target_time   = toc(eyeTrack_clock);
        
          Screen('FillRect', eyeTrack_window, grey, [0 0 screenXpixels screenYpixels], 0);  
          Screen('Flip', eyeTrack_window); 
        
          Screen('FillRect', stimulus_window, grey, [0 0 screenXpixels screenYpixels], 0);  
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
      
        if keyIsDown && KbName(keyCode) == 'q';
          is_running = 0;      
          close all;
          sca;
          expt_info.trial_records = trial_records;
          if save_records
            save(saveSTR,'expt_info');
          end
          
          break;
        end    
      end     
   
      if endCounter != 0  
        for i = 1 : (endCounter - startCounter + 1)
          trCount(i) = counter(startCounter + i - 1);
        end
      end
   
      % record trial data
      trial_records(tr_ind).stimulus = stimulus_image;
      trial_records(tr_ind).trCount = trCount;
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
        expt_info.trial_records = trial_records;
        if save_records
          save(saveSTR,'expt_info');
        end        
        break;
      end  
    
    
      if tr_ind >= max_trs
        is_running = 0;
        sca;     
      end   
    
      Datapixx('SetDoutValues', 0);
      Datapixx('RegWrRd');      
    end
  
catch 
  lasterr
  sca;  
end
end

  # last line