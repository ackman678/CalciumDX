%getWaveISI
%James Ackman, 1/7/2011
%wave ISI------------------------------------------------------------------
tmpres=region.timeres;
d=region.waveonsets;
e=region.waveoffsets;

d1 = [d size(region.traces,2)];
e1 = [0 e-d];
ints=diff([0 d1]) - e1;
if d(1) ~= 1
    ints=ints(2:end);
end
if d(end) ~= size(region.traces,2)
    ints=ints(1:end-1);
end
ints=ints*tmpres;

% disp(num2str(ints))