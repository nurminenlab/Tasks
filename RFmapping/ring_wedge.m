function trial_records = ring_wedge(debug_on)
  % prepare PsychoToolbox
  addpath('/home/vpixx/Tasks/Functions/');
  AssertOpenGL;
  sca;
  close all;

  # use mouse instead of eye tracker
  mouse_track = 0;

  animal = 'Sansa';
  saveSTR = ['/home/vpixx/MonkeyRecords/TrialRecords/',animal,'/','RF-map-wedge-trial_records-',date,'.mat'];
  STR = ['/home/vpixx/MonkeyRecords/TrialRecords/',animal,'/','RF-map-wedge-trial_records-',date];
  save_append = 0;
  while exist(saveSTR,'file') == 2
    save_append = save_append + 1;
    saveSTR = [STR,'_v',num2str(save_append),'.mat'];
  end  

  % user defined parameters
  if strcmp(animal,'Tully-')
    scaler = 0.6;
    
  elseif strcmp(animal,'Sansa')
    view_distance = 47;
    pix_per_cm = 36.2;
    va_in_pix  = va2pix(view_distance,pix_per_cm);
    Trans_mx_shift = [30 -30]; # a manual offset to the translation matrix of the eye tracker calibration. DEF in pixels. 
    
    fixation_target_deg = 1.5;    
    trackWin_deg = 2.5; # diameter
    ring_width_deg = 1;    
    
    waitframes = 7;    
    fill_fixation = 1;
    wait_fixation        = 0.75;
    rewardConsume_period = 2;
    ms                   = 10;    
    max_trs              = 10000;    
    max_fixation_duration = 10;
    min_fixation_duration = 0.1;
    reward_scaler = 0.5;
    FR = 120;    
    trackMarkerColor = [255,0,0];
    gaze_position = nan*ones(2,FR*ceil((wait_fixation+max_fixation_duration)));
  else
    fprintf('Strange animal \n');
  end  
  
  expt_info.view_distance = view_distance;
  expt_info.pix_per_cm = pix_per_cm;
  expt_info.va_in_pix = va_in_pix;
  expt_info.fixation_target_deg = fixation_target_deg;
  expt_info.ring_width_deg = ring_width_deg;
  expt_info.trackWin_deg = trackWin_deg;
  expt_info.waitframes = waitframes;
  expt_info.fill_fixation = fill_fixation;
  expt_info.wait_fixation = wait_fixation;
  expt_info.rewardConsume_period = rewardConsume_period;
  expt_info.max_trs  = max_trs;
  expt_info.max_fixation_duration = max_fixation_duration;  
  expt_info.min_fixation_duration = min_fixation_duration;
  expt_info.reward_scaler = reward_scaler;
  
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
  
  % Get the size of the on screen window, these are the same for both screens
  [screenXpixels, screenYpixels] = Screen('WindowSize',stimulus_window);
  
  %% NEW
  ifi1 = Screen('GetFlipInterval', eyeTrack_window);
  ifi2 = Screen('GetFlipInterval', stimulus_window);

  square_size = screenYpixels/2;
  if(fill_fixation == 1)
    [track_centerX, track_centerY] = WindowCenter(eyeTrack_window);
    fillStimulus = [track_centerX-square_size, track_centerY-square_size, track_centerX+square_size, track_centerY+square_size];
  else
    fillStimulus = [0 0 0 0];
  end
    
  wedge = randperm(180,90) + 180;
  wedgeCount = 1;
  ringCount = 1;
  
  ring = 1:1:16;
  ring = (va_in_pix*ring)/2; # /2 due to downstream    
    
  % Load marmoset face
  stimulus_image = 'face10.jpg';
  theImage = imread(stimulus_image);
  [s1, s2, s3] = size(theImage);
  
  % scale image rectangle
  rect = [1 1 fixation_target_deg*va_in_pix fixation_target_deg*va_in_pix];  
  eyePos_rect = [0 0 5 5];
  trackWindow_rect = [1 1 trackWin_deg*va_in_pix trackWin_deg*va_in_pix];
  trackWindow = (trackWin_deg*va_in_pix)/2;
  
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
  
  # grating
  windowPointer = stimulus_window;  
  
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
    
    
    if tr_ind > 1
      Screen('Flip', eyeTrack_window,0);
    end   
    
    % Grey screen (stay grey after exiting program)
    Screen('FillRect', stimulus_window, grey, rects(:,:,pos));
    Screen('FillRect', eyeTrack_window, grey, rects(:,:,pos));
    
    greyScreen_stimulus_vbl = Screen('Flip', stimulus_window);  
    greyScreen_eyeTrack_vbl = Screen('Flip', eyeTrack_window);  
    
    reward_consume_clock = tic();
    
    while toc(reward_consume_clock) < rewardConsume_period;
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
        expt_info.trial_records = trial_records;        
        save(saveSTR,'expt_info');        
        break;
      end
      
    end  
    
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
        # set digital out to 1 for channel XYZ
        Datapixx('SetDoutValues', 20);
        Datapixx('RegWrRd');
        on_target = 1;
        
        reaction_time = toc(wait_fixation_clock);
        eyeTrack_clock = tic();      
               
        trWedge = [];
        trRing  = []; 
        if(wedgeCount == 90)
          wedgeCount = 1;
          wedge = randperm(180, 90) + 180;
        else      
          wedgeCount = wedgeCount + 1;
        end
        
        inRing = [track_centerX - ring(ringCount), track_centerY - ring(ringCount), track_centerX + ring(ringCount), track_centerY + ring(ringCount)]; 
        outRing = [track_centerX - (ring(ringCount)+ring_width_deg*va_in_pix), ...
                    track_centerY - (ring(ringCount)+ring_width_deg*va_in_pix), ...
                    track_centerX + (ring(ringCount)+ring_width_deg*va_in_pix), ...
                    track_centerY + (ring(ringCount)+ring_width_deg*va_in_pix)];
        ringCount = ringCount + 1;
        
        if ringCount > length(ring)
          ringCount = 1;
          ring_ind = randperm(length(ring), length(ring));
          ring = ring(ring_ind);
        end
    
        turn = 0;
        
        if(mod(turn, 2) == 0)
          Screen('FillArc', eyeTrack_window, [0 0 0], fillStimulus, wedge(wedgeCount), 9);
          Screen('FillArc', stimulus_window, [0 0 0], fillStimulus, wedge(wedgeCount), 9);
        else
          Screen('FillOval', eyeTrack_window, [0 0 0], outRing);
          Screen('FillOval', eyeTrack_window, [128 128 128], inRing);
          Screen('FillOval', stimulus_window, [0 0 0], outRing);
          Screen('FillOval', stimulus_window, [128 128 128], inRing);
        end
        
        Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));
        Screen('FrameOval',eyeTrack_window,[0 0 255],trackWindow_rect(:,:,pos),3,3); # to mark fixation window
          
        Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos)); 
        Screen('FillOval', stimulus_window, [128 128 128], [0 0 0 0]);
        
        turn = turn + 1;
        trWedge(turn) = wedge(wedgeCount);
        trRing(turn) = ring(ringCount);
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
      
      reaction_time = NaN;
      trial_error = 'no_fixation';
      trWedge = NaN;
      trRing = NaN;
      fixation_duration = NaN;

      if keyIsDown && KbName(keyCode) == 'q';
        is_running = 0;
        Datapixx('StopAllSchedules');
        Datapixx('Close');      
        close all;
        sca;
        expt_info.trial_records = trial_records;        
        save(saveSTR,'expt_info');        
        break;
      end      
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
      
      if(wedgeCount == 90)
        wedgeCount = 1;
        wedge = randperm(180, 90) + 180;
      else      
        wedgeCount = wedgeCount + 1;
      end
      
      inRing = [track_centerX - ring(ringCount), track_centerY - ring(ringCount), track_centerX + ring(ringCount), track_centerY + ring(ringCount)]; 
      outRing = [track_centerX - (ring(ringCount)+ring_width_deg*va_in_pix), ...
                  track_centerY - (ring(ringCount)+ring_width_deg*va_in_pix), ...
                  track_centerX + (ring(ringCount)+ring_width_deg*va_in_pix), ...
                  track_centerY + (ring(ringCount)+ring_width_deg*va_in_pix)];
      ringCount = ringCount + 1;
        
      if ringCount > length(ring)
        ringCount = 1;
        ring_ind = randperm(length(ring), length(ring));
        ring = ring(ring_ind);
      end              
      
      
      if(mod(turn, 2) == 0)
        Screen('FillArc', eyeTrack_window, [0 0 0], fillStimulus, wedge(wedgeCount), 9);
        Screen('FillArc', stimulus_window, [0 0 0], fillStimulus, wedge(wedgeCount), 9);
      else
        Screen('FillOval', eyeTrack_window, [0 0 0], outRing);
        Screen('FillOval', eyeTrack_window, [128 128 128], inRing);
        Screen('FillOval', stimulus_window, [0 0 0], outRing);
        Screen('FillOval', stimulus_window, [128 128 128], inRing);
      end
      
      Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rects(:,:,pos));      
      Screen('FrameOval',eyeTrack_window,[0 0 255],trackWindow_rect(:,:,pos),3,3); # to mark fixation window
      
      Screen('DrawTexture', stimulus_window, stimulus_imageTexture, [], rects(:,:,pos));
      
      turn = turn + 1;
      trWedge(turn) = wedge(wedgeCount);
      trRing(turn)  = ring(ringCount);
      
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
      
      if sqrt((XY(1) - X(pos))^2 + (XY(2) - Y(pos))^2) >= trackWindow      
        # set digital out to 0 for channel XYZ
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        on_target = 0;        
        
        Screen('FillRect', eyeTrack_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', eyeTrack_window); 
        
        Screen('FillRect', stimulus_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', stimulus_window); 
        
        fixation_duration = toc(eyeTrack_clock);
        trial_error = 'broke_early';
        if fixation_duration > min_fixation_duration
          reward_size_time = reward_scaler*sqrt(fixation_duration);
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
      if toc(eyeTrack_clock) > max_fixation_duration
        # set digital out to 0 for channel XYZ
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        on_target = 0;
        fixation_duration = toc(eyeTrack_clock);
        
        Screen('FillRect', eyeTrack_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', eyeTrack_window); 
        
        Screen('FillRect', stimulus_window, grey, [0 2000 0 2000], 0);  
        Screen('Flip', stimulus_window); 
        
        fixation_duration = toc(eyeTrack_clock);# set digital out to 0 for channel XYZ
        reward_size_time = reward_scaler*sqrt(fixation_duration);
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
        expt_info.trial_records = trial_records;        
        save(saveSTR,'expt_info');        
        break;
      end    
    end  %   
    
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if keyIsDown && KbName(keyCode) == 'q';
      is_running = 0;      
      close all;
      sca; 
      expt_info.trial_records = trial_records;
      save(saveSTR,'expt_info');     
      break;
    end
    
    % record trial data
    trial_records(tr_ind).stimulus = stimulus_image;
    trial_records(tr_ind).trWedge = trWedge;
    trial_records(tr_ind).trRing = trRing;
    trial_records(tr_ind).positionX = X(pos);
    trial_records(tr_ind).positionY = Y(pos);
    trial_records(tr_ind).reaction_time  = reaction_time;
    trial_records(tr_ind).fixation_duration = fixation_duration;
    trial_records(tr_ind).trial_error = trial_error;
    trial_records(tr_ind).gaze_position = gaze_position; 
    
  end
