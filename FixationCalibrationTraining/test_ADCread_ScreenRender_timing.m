[FNAME, FPATH, FLTIDX] = uigetfile();
load([FPATH,FNAME]);
[bx,by] = compute_calibration_matrix(tr);

Scale_mx = eye(2);
Scale_mx(1) = bx(2);
Scale_mx(4) = by(2);
Trans_mx = [bx(1), by(1)]';

Datapixx('Open');

adcRate = 1e3;
nAdcSamples = 1000;

stimulus_duration = 1;
minreadFrames = 15;

[eyeTrack_window, eyeTrack_windowRect] = PsychImaging('OpenWindow', 0, 155);
Screen('BlendFunction', eyeTrack_window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

Datapixx('SetAdcSchedule', 0, adcRate, nAdcSamples, [0 2]);
Datapixx('StartAdcSchedule')
Datapixx('RegWrRd');

nreadFrames = 0;
tic();
while nreadFrames < minreadFrames
  Datapixx('RegWrRd');
  status = Datapixx('GetAdcStatus');
  nreadFrames = status.newBufferFrames;
  if (nreadFrames < minreadFrames)
    continue;
  else
    dada = Datapixx('ReadAdcBuffer', nreadFrames, -1);  
    # scale and translate
    dada = Scale_mx*dada + Trans_mx;
    Screen('DrawLines', eyeTrack_window, dada);
    t = toc();
    eyeTrack_vbl = Screen('Flip', eyeTrack_window);    
  end
end
KbStrokeWait(); 
sca();
