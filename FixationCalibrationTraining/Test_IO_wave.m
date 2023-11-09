AssertOpenGL;   % We use PTB-3
##client = tcp("172.27.85.73",5000);
##tcp_read(client,client.bytesavailable); # empty the buffer
##tcp_write(client,"OK"); # imager, we are ready!

Datapixx('Open');
Datapixx('StopAllSchedules');
Datapixx('RegWrRd');    % Synchronize DATAPixx registers to local register cache

% We'll make sure that all the TTL digital outputs are low before we start
Datapixx('SetDoutValues', 0);
Datapixx('RegWrRd');

% Define what we want a "trigger pulse" to look like,
% then download it to the DATAPixx.
% We'll arbitrarily say that it is 1 sample high, and 3 samples low.
doutWave = [16 0 0 0];
bufferAddress = 8e6;
Datapixx('WriteDoutBuffer', doutWave, bufferAddress);

% Define the schedule which will play the wave.
samplesPerTrigger = size(doutWave,2);
triggersPerFrame = 1;
samplesPerFrame = samplesPerTrigger * triggersPerFrame;
framesPerTrial = 600;       % We'll send triggers for 20 seconds at 120Hz
samplesPerTrial = samplesPerFrame * framesPerTrial;
#Datapixx('SetDoutSchedule', 0, [samplesPerFrame, 2], samplesPerTrial, bufferAddress, samplesPerTrigger);

% Insert your visual stimulus setup code here, finishing up with a Screen('Flip', window);
stimulus_window = Screen('OpenWindow',1,128);

% Tell the trigger schedule to start on the next vertical sync pulse
Datapixx('SetDoutSchedule', 0, [samplesPerFrame, 2], samplesPerTrial, bufferAddress, samplesPerTrigger);
for i = 1:200  
  Datapixx('StartDoutSchedule');
  Datapixx('RegWrVideoSync');
  Screen('Flip',stimulus_window,0);
  
  blank_clock = tic();
  while toc(blank_clock) < 11
    # blank screen
  endwhile
  disp(['Completed block number ', num2str(i)]);
  
  Datapixx('SetDoutSchedule', 0, [samplesPerFrame, 2], samplesPerTrial, bufferAddress, samplesPerTrigger);
end
