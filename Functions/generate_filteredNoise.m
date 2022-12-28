function [grText] = generate_filteredNoise(gridSize,windowPointer,parameters)
  
addpath('/home/vpixx/Tasks/Functions/');
echo off

SFfilter = Bandpass2(gridSize,parameters.lowCut_S,parameters.highCut_S);
ORfilter = OrientationBandpass(gridSize,parameters.orientation_low_S,parameters.orientation_high_S);

filter   = ORfilter.*SFfilter;
noise_patch = 2.*rand(gridSize,gridSize) -1;
ftd=filter.*fftshift(fft2(noise_patch));
iftd_noise = real(ifft2(ifftshift(ftd)));
iftd_noise = iftd_noise./max(abs(min(min(iftd_noise))),max(max(iftd_noise)));
iftd_noise = iftd_noise.*parameters.contrast_S;
raisedCosineMask = raised_cosine(parameters.plateauPixels_S,parameters.edgePixels_S,gridSize,gridSize,0,0,'R');

black = BlackIndex(windowPointer);  
white = WhiteIndex(windowPointer);  
gray = (black + white) / 2; 
if round(gray)==white
  gray=black;
end

absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);
grayScaleImageMatrix = gray + absoluteDifferenceBetweenWhiteAndGray*iftd_noise.*raisedCosineMask;
grText = Screen('MakeTexture', windowPointer, grayScaleImageMatrix);