# debug RFmap

addpath('/home/vpixx/Tasks/Functions/');
AssertOpenGL;
close all;
sca;

#PsychDebugWindowConfiguration;

scaler = 0.1;

if true
  Datapixx('Open')
  Datapixx('StopAllSchedules')
  Datapixx('DisableDacAdcLoopback');
  Datapixx('RegWrRd');
  adcRate = 1e3;
  baseBuffAddr = 5e6;
  minStreamFrames = 15;
  
  Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
  Datapixx('StartAdcSchedule');
  Datapixx('RegWrRd');
  # to force correct baseBuffAddr
  XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
  Datapixx('RegWrRd');  
end


screens = Screen('Screens');
eyeTrack_screenNumber = min(screens);
stimulus_screenNumber = max(screens);

white = WhiteIndex(stimulus_screenNumber);
black = BlackIndex(stimulus_screenNumber);
grey = 128;
inc = white - grey;

eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,grey);
[screenXpixels, screenYpixels] = Screen('WindowSize',eyeTrack_window);

stimulus_image = 'face10.jpg';
theImage = imread(stimulus_image);
[s1, s2, s3] = size(theImage);

rect = [0 0 s1*scaler s2*scaler];
rect = CenterRectOnPoint(rect, 960,540);
eyeTrack_imageTexture = Screen('MakeTexture', eyeTrack_window, theImage);
Screen('DrawTexture', eyeTrack_window, eyeTrack_imageTexture, [], rect);
Screen('Flip', eyeTrack_window, 0,1);
KbStrokeWait();
sca;
close all;
