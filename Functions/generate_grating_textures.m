function grText = generate_grating_textures(gridSize,orientations,pixelsPerPeriod,plateauPixels,edgePixels,windowPointer,screenNumber,contrast)
  
addpath('/home/vpixx/Tasks/Functions/');
echo off

orientations = orientations * pi / 180; 
spatialFrequency = 1 / pixelsPerPeriod; 
radiansPerPixel = spatialFrequency * (2 * pi);

raisedCosineMask = raised_cosine(plateauPixels,edgePixels,gridSize,gridSize,0,0,'R');

black = BlackIndex(windowPointer);  
white = WhiteIndex(windowPointer);  
gray = (black + white) / 2;  
if round(gray)==white
  gray=black;
end

absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);
[x y] = meshgrid(linspace(-gridSize/2,gridSize/2,gridSize), linspace(-gridSize/2,gridSize/2,gridSize));

for i = 1:length(orientations)
  a=cos(orientations(i))*radiansPerPixel;
  b=sin(orientations(i))*radiansPerPixel;  
  gratingMatrix(:,:,i) = contrast*sin(a*x+b*y);
  grayScaleImageMatrix(:,:,i) = gray + absoluteDifferenceBetweenWhiteAndGray*gratingMatrix(:,:,i).*raisedCosineMask;
end

grText = NaN * ones(size(orientations));
for i = 1:size(grayScaleImageMatrix,3)
  grText(i) = Screen('MakeTexture', windowPointer, round(grayScaleImageMatrix(:,:,i)));
end

