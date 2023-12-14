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
doutWave = zeros(1,10);
doutWave(1) = 16;
doutWave = repmat(doutWave,1,600);
samplesPerTrial = size(doutWave,2);
bufferAddress = 8e6;
Datapixx('WriteDoutBuffer', doutWave, bufferAddress);

% Insert your visual stimulus setup code here, finishing up with a Screen('Flip', window);
stimulus_window = Screen('OpenWindow',1,128);

% Tell the trigger schedule to start on the next vertical sync pulse
Datapixx('SetDoutSchedule', 0, 1000, samplesPerTrial, bufferAddress);
for i = 1:1
  Datapixx('StartDoutSchedule');
  Datapixx('RegWrVideoSync');
  Screen('Flip',stimulus_window,0);
  status = Datapixx('GetDoutStatus');
  while status.scheduleRunning == 1
    # wait the schedule to finish
  endwhile
  blank_clock = tic();
  while toc(blank_clock) < 5
    # blank screen
  endwhile
  Datapixx('SetDoutSchedule', 0, 1000, samplesPerTrial, bufferAddress);
end
