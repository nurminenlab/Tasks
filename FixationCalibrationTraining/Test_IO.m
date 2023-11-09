function testing = Test_IO(debug_on)
  Datapixx('Open');
  adcRate = 1e3;
  baseBuffAddr = 0;
  Datapixx('SetAdcSchedule',0,adcRate,0,[0 2],baseBuffAddr,adcRate);
  Datapixx('StartAdcSchedule');
  Datapixx('RegWrRd');
  # to force correct baseBuffAddr
  XY = Datapixx('ReadAdcBuffer', 1, baseBuffAddr);
  Datapixx('RegWrRd');

  Datapixx('SetDoutValues', 16);
  Datapixx('RegWrRd');
  
  a = tic();      
  while toc(a) < 0.3
  end  
  #if toc(a) > 0.5
            
  #end
      
  Datapixx('SetDoutValues', 0);
  Datapixx('RegWrRd');
  #Datapixx('Close');
end