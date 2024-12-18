AssertOpenGL;   % We use PTB-3

% Open Datapixx, and stop any schedules which might already be running
Datapixx('Open');
Datapixx('StopAllSchedules');
Datapixx('RegWrRd');    % Synchronize DATAPixx registers to local register cache

% Show the output voltage range for each DAC channel
dacRanges = Datapixx('GetDacRanges');

KbStrokeWait();
sca;
