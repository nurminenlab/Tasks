function [grText] = generate_filteredNoise(gridSize,windowPointer,parameters)
  
addpath('/home/vpixx/Tasks/Functions/');
echo off

SFfilter_surround = Bandpass2(gridSize,parameters.lowCut_S,parameters.highCut_S);
ORfilter_surround = OrientationBandpass(gridSize,parameters.orientation_low_S,parameters.orientation_high_S);
SFfilter_center = Bandpass2(gridSize,parameters.lowCut_C,parameters.highCut_C);
ORfilter_center = OrientationBandpass(gridSize,parameters.orientation_low_C,parameters.orientation_high_C);

raisedCosineMask = raised_cosine(parameters.plateauPixels_S,parameters.edgePixels_S,gridSize,gridSize,0,0,'R');
raisedCosineMask_center  = raised_cosine(parameters.plateauPixels_C,parameters.edgePixels_C,gridSize,gridSize,0,0,'R');
raisedCosineMask_annulus = raised_cosine(parameters.plateauPixels_C,parameters.edgePixels_C,gridSize,gridSize,0,0,'H');

surround_filter = ORfilter_surround.*SFfilter_surround;
surround_noise_patch = 2.*rand(gridSize,gridSize) -1;
surround_ftd = surround_filter.*fftshift(fft2(surround_noise_patch));
surround_noise = real(ifft2(ifftshift(surround_ftd)));
surround_noise = surround_noise./max(abs(min(min(surround_noise))),max(max(surround_noise)));
surround_noise = surround_noise.*parameters.contrast_S;
surround_noise = surround_noise.*raisedCosineMask.*raisedCosineMask_annulus;

center_filter = ORfilter_center.*SFfilter_center;
center_noise_patch = 2.*rand(gridSize,gridSize) -1;
center_ftd = center_filter.*fftshift(fft2(center_noise_patch));
center_noise = real(ifft2(ifftshift(center_ftd)));
center_noise = center_noise./max(abs(min(min(center_noise))),max(max(center_noise)));
center_noise = center_noise.*parameters.contrast_C;
center_noise = center_noise.*raisedCosineMask_center;
grayScaleImageMatrix = center_noise + surround_noise;

black = BlackIndex(windowPointer);  
white = WhiteIndex(windowPointer);  
gray = (black + white) / 2; 
if round(gray)==white
  gray=black;
end

absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);
grayScaleImageMatrix = gray + absoluteDifferenceBetweenWhiteAndGray.*grayScaleImageMatrix;
grText = Screen('MakeTexture', windowPointer, grayScaleImageMatrix);