function GratingDemo
addpath('/home/vpixx/Tasks/Functions/');
% GratingDemo
%
% Displays a stationary grating.  See also DriftDemo, DriftDemo2, DriftDemo3 and DriftWaitDemo.

% ---------- Program History ----------

% 07/01/1999 dgp Added arbitrary orientation.
% 12/10/2001 awi Added font conditional.
% 02/21/2002 dgp Mentioned DriftDemo.
% 04/03/2002 awi Merged OS9 and Win versions, which had fallen out of sync. 
% 04/13/2002 dgp Used Arial, eliminating need for conditional.
% 07/15/2003 dgp Added comments explaining f and lambda.
% 08/16/2006 rhh Added user-friendly parameters, such as tiltInDegrees,
%                pixelsPerPeriod, periodsCoveredByOneStandardDeviation and widthOfGrid.
% 08/18/2006 rhh Expanded comments and created comment sections.
% 10/04/2006 dhb Minimize warnings.
% 10/11/2006 dhb Use maximum available screen.
% 10/14/2006 dhb Save and restore altered prefs, more extensive comments for them
% 07/12/2006 prf Changed method of rotating the grating

% ---------- Parameter Setup ----------
% Initializes the program's parameters.

% Prevents MATLAB from reprinting the source code when the program runs.
echo off
stimulus_image = 'face8.jpg';
theImage = imread(stimulus_image);

% *** To rotate the grating, set tiltInDegrees to a new value.
tiltInDegrees = [0,90]; % The tilt of the grating in degrees.
tiltInRadians = tiltInDegrees * pi / 180; % The tilt of the grating in radians.

% *** To lengthen the period of the grating, increase pixelsPerPeriod.
pixelsPerPeriod = 33; % How many pixels will each period/cycle occupy?
spatialFrequency = 1 / pixelsPerPeriod; % How many periods/cycles are there in a pixel?
radiansPerPixel = spatialFrequency * (2 * pi); % = (periods per pixel) * (2 pi radians per period)

plateauCycles = 3;
plateauPixels = plateauCycles*pixelsPerPeriod;

edgeCycles = 0.4;
edgePixels = edgeCycles*pixelsPerPeriod;

% *** If the grating is clipped on the sides, increase widthOfGrid.
widthOfGrid = 512;
halfWidthOfGrid = widthOfGrid / 2;
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  % widthArray is used in creating the meshgrid.

raisedCosineMask = raised_cosine(plateauPixels,edgePixels,widthOfGrid,widthOfGrid,0,0,'R');
% For an explanation of the try-catch block, see the section "Error Handling"
% at the end of this document.
try
  % ---------- Window Setup ----------
  % Opens a window.

  % Screen is able to do a lot of configuration and performance checks on
  % open, and will print out a fair amount of detailed information when
  % it does.  These commands supress that checking behavior and just let
  % the demo go straight into action.  See ScreenTest for an example of
  % how to do detailed checking.
  oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
  oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);

  % Find out how many screens and use largest screen number.
  whichScreen = max(Screen('Screens'));

  % Opens a graphics window on the main monitor (screen 0).  If you have
  % multiple monitors connected to your computer, then you can specify
  % a different monitor by supplying a different number in the second
  % argument to OpenWindow, e.g. Screen('OpenWindow', 2).
  window = Screen('OpenWindow', whichScreen);

  % Hides the mouse cursor
  HideCursor;

  % ---------- Color Setup ----------
  % Gets color values.

  % Retrieves color codes for black and white and gray.
  black = BlackIndex(window);  % Retrieves the CLUT color code for black.
  white = WhiteIndex(window);  % Retrieves the CLUT color code for white.
  gray = (black + white) / 2;  % Computes the CLUT color code for gray.
  if round(gray)==white
    gray=black;
  end

  % Taking the absolute value of the difference between white and gray will
  % help keep the grating consistent regardless of whether the CLUT color
  % code for white is less or greater than the CLUT color code for black.
  absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

  % ---------- Image Setup ----------
  % Stores the image in a two dimensional matrix.

  % Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
  % the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
  % and y = y(x0, y0) corresponds to the y-coordinate of element "i"
  [x y] = meshgrid(linspace(-widthOfGrid/2,widthOfGrid/2,widthOfGrid), linspace(-widthOfGrid/2,widthOfGrid/2,widthOfGrid));
 
  % Replaced original method of changing the orientation of the grating
  % (gradient = y - tan(tiltInRadians) .* x) with sine and cosine (adapted from DriftDemo). 
  % Use of tangent was breakable because it is undefined for theta near pi/2 and the period
  % of the grating changed with change in theta.  

  for i = 1:length(tiltInRadians)
    a=cos(tiltInRadians(i))*radiansPerPixel;
    b=sin(tiltInRadians(i))*radiansPerPixel;

  % Converts meshgrid into a sinusoidal grating, where elements

  % along a line with angle theta have the same value and where the
  % period of the sinusoid is equal to "pixelsPerPeriod" pixels.
  % Note that each entry of gratingMatrix varies between minus one and
  % one; -1 <= gratingMatrix(x0, y0)  <= 1
    gratingMatrix(:,:,i) = sin(a*x+b*y);
    grayScaleImageMatrix(:,:,i) = gray + absoluteDifferenceBetweenWhiteAndGray*gratingMatrix(:,:,i).*raisedCosineMask;
  end  
  
  % ---------- Image Display ---------- 
  % Displays the image in the window.

  % Colors the entire window gray.
  Screen('FillRect', window, gray);

  % Writes the image to the window.
  grText = NaN * ones(size(grayScaleImageMatrix,3));
  for i = 1:size(grayScaleImageMatrix,3)
    grText(i) = Screen('MakeTexture', window, grayScaleImageMatrix(:,:,i));
  end  
 
  imPointer = Screen('MakeTexture', window, theImage);
 
  Screen('DrawTexture', window, imPointer);
  Screen('DrawTexture', window, grText(2));
  
  
  % Writes text to the window.
  currentTextRow = 0;
  Screen('DrawText', window, sprintf('black = %d, white = %d', black, white), 0, currentTextRow, black);
  currentTextRow = currentTextRow + 20;
  Screen('DrawText', window, 'Press any key to exit.', 0, currentTextRow, black);

  % Updates the screen to reflect our changes to the window.
  Screen('Flip', window);

  % Waits for the user to press a key.
  KbWait;

  % ---------- Window Cleanup ---------- 

  % Closes all windows.
  sca;

  % Restores the mouse cursor.
  ShowCursor;

  % Restore preferences
  Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
  Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
catch
  % ---------- Error Handling ---------- 
  % If there is an error in our code, we will end up here.

  % The try-catch block ensures that Screen will restore the display and return us
  % to the MATLAB prompt even if there is an error in our code.  Without this try-catch
  % block, Screen could still have control of the display when MATLAB throws an error, in
  % which case the user will not see the MATLAB prompt.
  sca;

  % Restores the mouse cursor.
  ShowCursor;

  % Restore preferences
  Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
  Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

  % We throw the error again so the user sees the error description.
  psychrethrow(psychlasterror);
end
