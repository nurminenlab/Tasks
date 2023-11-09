% prepare PsychoToolbox
AssertOpenGL;
sca;
close all;

#PsychDebugWindowConfiguration;

# general parameters
mouse_track = 0;
always_centered_mouse = false;
always_centered_tracker = false;
max_trials  = 100;
trackWindow = 60;

Datapixx('Open');
adcRate = 1e3;
baseBuffAddr = 0;
minStreamFrames = 15;

framesPerTrial = 660;

% stimulus parameters
f=0.02;
cyclespersecond=1;
angle=0;
angle2=90;
movieDurationSecs=5; 
texsize=1920; 

stimulus_image = 'face10.jpg';

# no adjustable parameters beyond this line
if ~mouse_track 
  [FNAME, FPATH, FLTIDX] = uigetfile();
  load([FPATH,FNAME]);
  [bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr,1);
  Scale_mx = eye(2);
  Scale_mx(1) = bx(2);
  Scale_mx(4) = by(2);
  Trans_mx = [bx(1), by(1)]';      
end

% initialize windows, draw stimuli etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doutWave = [16 0 0 0];
bufferAddress = 8e6;
Datapixx('WriteDoutBuffer', doutWave, bufferAddress);
Datapixx('RegWr');
reward_bufferAddress = 5e6;

samplesPerTrigger = size(doutWave,2);
triggersPerFrame = 1;
samplesPerFrame = samplesPerTrigger * triggersPerFrame;
samplesPerTrial = samplesPerFrame * framesPerTrial;

theImage = imread(stimulus_image);

screenNumber   = max(Screen('Screens'));
iTscreenNumber = min(Screen('Screens'));

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);

if gray == white
  gray=white / 2;
end
inc=white-gray;

w   = Screen('OpenWindow',screenNumber, gray);
iTw = Screen('OpenWindow',iTscreenNumber, gray);
[X,Y] = WindowCenter(w);
AssertGLSL;

stimulus_imageTexture = Screen('MakeTexture', w, theImage);
iT_imageTexture       = Screen('MakeTexture', iTw, theImage);

Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('BlendFunction', iTw, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

glsl   = MakeTextureDrawShader(w, 'SeparateAlphaChannel');
iTglsl = MakeTextureDrawShader(iTw, 'SeparateAlphaChannel');

p=ceil(1/f); 
fr=f*2*pi;
visiblesize=2*texsize+1;
x = meshgrid(-texsize:texsize + p, -texsize:texsize);
grating = gray + inc*cos(fr*x);
[x,y]=meshgrid(-texsize:texsize, -texsize:texsize);
circle = white * (x.^2 + y.^2 <= (texsize)^2);
grating(:,:,2) = white;
gratingtex1  = Screen('MakeTexture', w, grating , [], [], [], [], glsl); 
iTgratingtex1 = Screen('MakeTexture', iTw, grating , [], [], [], [], iTglsl);

srcRect=[0 0 visiblesize visiblesize];    
ifi=Screen('GetFlipInterval', w);
waitframes = 1;
waitduration = waitframes * ifi;    
p = 1/f;
shiftperframe = cyclespersecond * p * waitduration;

trackMarkerColor = [255,0,0];
eyePos_rect = [0 0 5 5];
trackWindow_rect = [1 1 trackWindow*2 trackWindow*2];
trackWindow_rect = CenterRectOnPoint(trackWindow_rect,X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mouse_track
  Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
  Datapixx('StartAdcSchedule');
  Datapixx('RegWrRd');
  # to force correct baseBuffAddr
  XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
  Datapixx('RegWrRd');
else
  XY = ones(2,1)*nan;
end

# open Datapixx
Datapixx('Open')
adcRate = 1e3;
baseBuffAddr = 0;
minStreamFrames = 15;

is_running = true;
on_target  = false; 
tr_ind     = 0;  
tracking   = false;
HideCursor(iTw);
while is_running

  tr_ind = tr_ind + 1;
  if rem(tr_ind,2) == 0
    mplier = -1;
  else
    mplier = 1;
  endif
  vbl = Screen('Flip', w);  

  # check keyboard
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if KbName(keyCode) == 'q'
    is_running = false;
    sca;
    break;    
  endif

  ######################
  # blank loop 
  vblendtime = vbl + movieDurationSecs;
  blank_clock = tic();
  while (vbl < vblendtime) && is_running
    
    # check keyboard
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if KbName(keyCode) == 'q'
      is_running = false;
      sca;
      break;    
    endif
    
    # TO DO asynchronous tracking
    if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else        
          if always_centered_tracker
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2); # for consistency
            [x,y] = WindowCenter(iTw);
            XY(1) = x;
            XY(2) = y;
          else            
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
            XY = Scale_mx*XY + Trans_mx;
          end          
        end
    else
        if always_centered_mouse
          [x,y] = WindowCenter(iTw);
        else          
          [x,y,buttons,focus] = GetMouse(iTw);
        endif        
        XY(1) = x;
        XY(2) = y;
    end
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    
    if on_target || (mod(vbl,2) < 1)
      #Screen('DrawTexture', w, stimulus_imageTexture);
      #Screen('DrawTexture', iTw, iT_imageTexture);            
    end;
    
    Screen('FillOval', iTw, trackMarkerColor, eyePos_rect);
    Screen('FrameOval',iTw, [0 0 255], trackWindow_rect, 3,3);
    if toc(blank_clock) >= (movieDurationSecs - 0.5)
      Datapixx('SetDoutSchedule', 0, [samplesPerFrame, 2], samplesPerTrial, bufferAddress, samplesPerTrigger);
      Datapixx('StartDoutSchedule');
      Datapixx('RegWrVideoSync');
    endif
    
    Screen('Flip', iTw,vbl + (waitframes - 0.5) * ifi);
    vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);

    if sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) < trackWindow
      on_target = true;
      if ~tracking 
        fix_clock = tic();
        tracking = true;
      endif
    end

    if tracking && (sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) >= trackWindow)
      disp('reward');
      on_target = false;
      tracking  = false;
      fix_time  = toc(fix_clock);
##      reward_size_time = 0.25*sqrt(fix_time);
##      reward_samples = floor(reward_size_time*1000);
##      reward_doutWave = [ones(1,reward_samples),0];
##      Datapixx('WriteDoutBuffer', reward_doutWave, reward_bufferAddress);
##      Datapixx('SetDoutSchedule', 0, 1000, reward_samples+1, reward_bufferAddress,reward_samples+1);
##      Datapixx('RegWr');      
    end        
  end;
  
  
  # vertical grating loop
  i=0;
  vblendtime = vbl + movieDurationSecs;
  while (vbl < vblendtime) && is_running
    
    # check keyboard
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if KbName(keyCode) == 'q'
      is_running = false;
      sca;
      break;    
    endif
    
    # TO DO asynchronous tracking
    if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else
          if always_centered_tracker
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2); # for consistency
            [x,y] = WindowCenter(iTw);
            XY(1) = x;
            XY(2) = y;
          else            
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
            XY = Scale_mx*XY + Trans_mx;
         end
        end       
    else
        if always_centered_mouse
          [x,y] = WindowCenter(iTw);
        else
          [x,y,buttons,focus] = GetMouse(iTw);  
        endif        
        XY(1) = x;
        XY(2) = y;
    end
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));    
    
    yoffset = mod(i*shiftperframe,p);
    i=i+1;
    
    Screen('DrawTexture', w, gratingtex1, srcRect, [], angle, [], [], [], [], [], [0, mplier*yoffset, 0, 0]);       
    Screen('DrawTexture', iTw, iTgratingtex1, srcRect, [], angle, [], [], [], [], [], [0, mplier*yoffset, 0, 0]);        
    
    if on_target || (mod(vbl,2) < 1)
      #Screen('DrawTexture', w, stimulus_imageTexture);
      #Screen('DrawTexture', iTw, iT_imageTexture);
    end;        
    
    Screen('FillOval', iTw, trackMarkerColor, eyePos_rect);
    Screen('FrameOval',iTw, [0 0 255], trackWindow_rect, 3,3);
    Screen('Flip', iTw,vbl + (waitframes - 0.5) * ifi);
    vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
    
    if sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) < trackWindow
      on_target = true;
      if ~tracking 
        fix_clock = tic();
        tracking = true;
      endif
    end

    if tracking && (sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) >= trackWindow)
      disp('reward');
      on_target = false;
      tracking  = false;
      fix_time  = toc(fix_clock);
      reward_size_time = 0.25*sqrt(fix_time);
##      reward_samples = floor(reward_size_time*1000);
##      reward_doutWave = [ones(1,reward_samples),0];
##      Datapixx('WriteDoutBuffer', reward_doutWave, reward_bufferAddress);
##      Datapixx('SetDoutSchedule', 0, 1000, reward_samples+1, reward_bufferAddress);
##      Datapixx('RegWr');      
    end
    
  end;
  
  # check keyboard
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if KbName(keyCode) == 'q'
    is_running = false;
    sca;    
    break;
  endif  
  
  # blank loop 
  vblendtime = vbl + movieDurationSecs;
  blank_clock = tic();
  while (vbl < vblendtime) && ~KbCheck    

    # check keyboard
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if KbName(keyCode) == 'q'
      is_running = false;
      sca;
      break;    
    endif
    
    # TO DO asynchronous tracking
    if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else  
          if always_centered_tracker
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2); # for consistency
            [x,y] = WindowCenter(iTw);
            XY(1) = x;
            XY(2) = y;
          else            
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);
            XY = Scale_mx*XY + Trans_mx;
         end
        end
      else
        if always_centered_mouse
          [x,y] = WindowCenter(iTw);
        else
          [x,y,buttons,focus] = GetMouse(iTw);  
        endif        
        XY(1) = x;
        XY(2) = y;
    end
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));
    
    if on_target || (mod(vbl,2) < 1)
##      Screen('DrawTexture', w, stimulus_imageTexture);
##      Screen('DrawTexture', iTw, iT_imageTexture);
    end;

    Screen('FillOval', iTw, trackMarkerColor, eyePos_rect);
    Screen('FrameOval',iTw, [0 0 255], trackWindow_rect, 3,3);
    
    if toc(blank_clock) >= (movieDurationSecs - 0.5)
      Datapixx('SetDoutSchedule', 0, [samplesPerFrame, 2], samplesPerTrial, bufferAddress, samplesPerTrigger);
      Datapixx('StartDoutSchedule');
      Datapixx('RegWrVideoSync');
    endif
    
    Screen('Flip', iTw,vbl + (waitframes - 0.5) * ifi);
    vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);

    if sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) < trackWindow
      on_target = true;
      if ~tracking 
        fix_clock = tic();
        tracking = true;
      endif
    end

    if tracking && (sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) >= trackWindow)
      disp('reward');
      on_target = false;
      tracking  = false;
      fix_time  = toc(fix_clock);
##      reward_size_time = 0.25*sqrt(fix_time);
##      reward_samples = floor(reward_size_time*1000);
##      reward_doutWave = [ones(1,reward_samples),0];
##      Datapixx('WriteDoutBuffer', reward_doutWave, reward_bufferAddress);
##      Datapixx('SetDoutSchedule', 0, 1000, reward_samples+1, reward_bufferAddress);
##      Datapixx('RegWr');      
    end
    
  end;     
  
  # check keyboard
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if KbName(keyCode) == 'q'
    is_running = false;
    sca;
    break;    
  endif
  
  # horizontal grating loop
  vblendtime = vbl + movieDurationSecs;
  while (vbl < vblendtime) && is_running
    
    # check keyboard
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if KbName(keyCode) == 'q'
      is_running = false;
      sca;
      break;    
    endif
    
    # TO DO asynchronous tracking
    if ~mouse_track
        Datapixx('RegWrRd');
        status = Datapixx('GetAdcStatus');
        nReadFrames = status.newBufferFrames;   
        if (nReadFrames < minStreamFrames)
          continue;
        else        
          if always_centered_tracker
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2); # for consistency
            [x,y] = WindowCenter(iTw);
            XY(1) = x;
            XY(2) = y;
          else            
            XY = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
            XY = Scale_mx*XY + Trans_mx;
         end
        end
      else
        if always_centered_mouse
          [x,y] = WindowCenter(iTw);
        else
          [x,y,buttons,focus] = GetMouse(iTw);
        endif        
        XY(1) = x;
        XY(2) = y;
    end
    eyePos_rect = CenterRectOnPoint(eyePos_rect,XY(1),XY(2));    
    
    yoffset = mod(i*shiftperframe,p);
    i=i+1;
    Screen('DrawTexture', w, gratingtex1, srcRect, [], angle2, [], [], [], [], [], [0, mplier*yoffset, 0, 0]);       
    Screen('DrawTexture', iTw, iTgratingtex1, srcRect, [], angle2, [], [], [], [], [], [0, mplier*yoffset, 0, 0]);
    
    if on_target || (mod(vbl,2) < 1)
##      Screen('DrawTexture', w, stimulus_imageTexture);
##      Screen('DrawTexture', iTw, iT_imageTexture);
    end;
    
    Screen('FillOval', iTw, trackMarkerColor, eyePos_rect);
    Screen('FrameOval',iTw, [0 0 255], trackWindow_rect, 3,3);    
    Screen('Flip', iTw,vbl + (waitframes - 0.5) * ifi);
    vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);

    if sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) < trackWindow
      on_target = true;
      if ~tracking 
        fix_clock = tic();
        tracking = true;
      endif
    end

    if tracking && (sqrt((XY(1) - X)^2 + (XY(2) - Y)^2) >= trackWindow)
      disp('reward');
      on_target = false;
      tracking  = false;
      fix_time  = toc(fix_clock);
##      reward_size_time = 0.25*sqrt(fix_time);
##      reward_samples = floor(reward_size_time*1000);
##      reward_doutWave = [ones(1,reward_samples),0];
##      Datapixx('WriteDoutBuffer', reward_doutWave, reward_bufferAddress);
##      Datapixx('SetDoutSchedule', 0, 1000, reward_samples+1, reward_bufferAddress);
##      Datapixx('RegWr');     
    end    
    
  end;
  
  # check keyboard
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  if KbName(keyCode) == 'q'
    is_running = false;
    sca;
    break;    
  endif  
  
  if tr_ind == max_trials
    is_running = false;
    sca;
  endif
  
end


