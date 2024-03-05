AssertOpenGL;
sca;
#PsychDebugWindowConfiguration;
% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if available
eyeTrack_screenNumber = min(screens);


% Define black and white
white = WhiteIndex(eyeTrack_screenNumber);
black = BlackIndex(eyeTrack_screenNumber);
grey = 128;
inc = white - grey;

% Open an on screen window
eyeTrack_window = Screen('OpenWindow',eyeTrack_screenNumber,grey);
Screen('GetFlipInfo',eyeTrack_window,1);
Screen('FillRect',eyeTrack_window,grey);
Screen('Flip',eyeTrack_window);
flipinfo = Screen('GetFlipInfo',eyeTrack_window);
sca;


