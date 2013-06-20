%getWaveDur
%James Ackman, 1/7/2011
%wave durations-- first onset to last onset----
newdur = region.timeres*(region.waveoffsets-region.waveonsets);
disp(num2str(newdur))