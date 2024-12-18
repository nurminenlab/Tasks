function trial_records = wedge_test(debug_on)
  % prepare PsychoToolbox
  addpath('/home/vpixx/Tasks/Functions/');
  AssertOpenGL;
  sca;
  close all;
  
  wedge_angle = 8.75;     
  wedge = 180 + 90;   

  if debug_on;
    PsychDebugWindowConfiguration;
  end
  
  % Get the screen numbers
  screens = Screen('Screens');
  
  % Draw to the external screen if available
  eyeTrack_screenNumber = min(screens);      
  grey = 128;  
  
  % Open an on screen window  
  eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,grey);
  
  % Get the size of the on screen window, these are the same for both screens
  [screenXpixels, screenYpixels] = Screen('WindowSize',eyeTrack_window);
  square_size = screenYpixels/2;
  [track_centerX, track_centerY] = WindowCenter(eyeTrack_window);
  fillStimulus = [track_centerX-square_size, track_centerY-square_size, track_centerX+square_size, track_centerY+square_size];
  
  Screen('FillArc', eyeTrack_window, [0 0 0], fillStimulus, wedge, wedge_angle);
  Screen('Flip', eyeTrack_window);  
  
  KbStrokeWait();
  sca;