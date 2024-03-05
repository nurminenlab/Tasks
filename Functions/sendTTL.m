function sendTTL(upvalue, downvalue, time_secs)
  # send TTL
  now = tic();
  Datapixx('SetDoutValues', upvalue);
  Datapixx('RegWrRd');
  while toc(now) < time_secs
    # do nothing
  end
  Datapixx('SetDoutValues', downvalue);          
  Datapixx('RegWrRd');  
endfunction
