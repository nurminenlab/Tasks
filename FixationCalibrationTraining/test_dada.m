% user defined parameters
scaler = 0.5;
stimulus_duration = 1;
isi_duration = 1;
ms = 10;

# load fixation calibration
[FNAME, FPATH, FLTIDX] = uigetfile();
load([FPATH,FNAME]);
[bx,by, HV, VV, HP, VP] = compute_calibration_matrix(tr);
Scale_mx = eye(2);
Scale_mx(1) = bx(2);
Scale_mx(4) = by(2);
Trans_mx = [bx(1), by(1)]';

% set-up Datapixx
Datapixx('Open');
adcRate = 1e3;
nAdcSamples = floor((stimulus_duration + isi_duration) * adcRate);
minStreamFrames = 15;

Datapixx('SetAdcSchedule', 0, adcRate, nAdcSamples, [0 2]);
Datapixx('StartAdcSchedule')
tic();
while toc() < isi_duration;
  [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
  % Start eye tracking and flip the stimulus to the screen        
  Datapixx('RegWrRd');
  status = Datapixx('GetAdcStatus');
  nReadFrames = status.newBufferFrames;   
  if (nReadFrames < minStreamFrames)
    continue;
   else
    dada = median(Datapixx('ReadAdcBuffer', nReadFrames, -1),2);         
    dada = Scale_mx*dada + Trans_mx;    
  end
end

Datapixx('StopAllSchedules');
Datapixx('Close');
  