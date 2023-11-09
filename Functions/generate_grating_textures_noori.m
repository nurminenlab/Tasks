function [grText,pixelsPerPeriodLUT,plateauPixelsLUT,contrastLUT] = generate_grating_textures_noori(gridSize,pixelsPerPeriod,plateauPixels,contrast,edgePixels,windowPointer,screenNumber,shader)
  
addpath('/home/vpixx/Tasks/Functions/');
echo off

orientation = 0;
black = BlackIndex(windowPointer);  
white = WhiteIndex(windowPointer);  
gray = (black + white) / 2;  
if round(gray)==white
  gray=black;
end
absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

[x y] = meshgrid(linspace(-gridSize/2,gridSize/2,gridSize), linspace(-gridSize/2,gridSize/2,gridSize));
grayScaleImageMatrix = ones(gridSize, gridSize,2);
txt = 0;
grText = NaN*ones(1,length(pixelsPerPeriod)*length(plateauPixels)*length(contrast));

pixelsPerPeriodLUT = NaN*ones(1,length(pixelsPerPeriod)*length(plateauPixels)*length(contrast));
plateauPixelsLUT   = NaN*ones(1,length(pixelsPerPeriod)*length(plateauPixels)*length(contrast));
contrastLUT        = NaN*ones(1,length(pixelsPerPeriod)*length(plateauPixels)*length(contrast));

for per = 1:length(pixelsPerPeriod)
  for plat = 1:length(plateauPixels)
    for c = 1:length(contrast)

      txt = txt + 1;
      spatialFrequency = 1 / pixelsPerPeriod(per); 
      radiansPerPixel = spatialFrequency * (2 * pi);
      raisedCosineMask = white * raised_cosine(plateauPixels(plat),edgePixels,gridSize,gridSize,0,0,'R');
      a=cos(orientation)*radiansPerPixel;
      b=sin(orientation)*radiansPerPixel;  
      gratingMatrix = contrast(c)*sin(a*x+b*y);

      # mask in the alpha channel      
      grayScaleImageMatrix(:,:,1) = round(gray + absoluteDifferenceBetweenWhiteAndGray*gratingMatrix);
      grayScaleImageMatrix(:,:,2) = raisedCosineMask;

      #grText = NaN * ones(size(orientations));
      grText(txt) = Screen('MakeTexture', windowPointer, grayScaleImageMatrix, [], [], [], [], shader);
      pixelsPerPeriodLUT(txt) = pixelsPerPeriod(per);
      plateauPixelsLUT(txt)   = plateauPixels(plat);
      contrastLUT(txt)        = contrast(c);
    end  
  end
end


