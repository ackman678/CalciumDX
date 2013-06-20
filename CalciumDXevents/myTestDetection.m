figure();
subplot(7,1,1)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxdettrial(trSign*region.traces(num,:)); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,2)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP(trSign*region.traces,num,'no'); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,3)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_FastAndSlow(trSign*region.traces,num,'no'); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,4)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHannFast(trSign*region.traces,num,'no'); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,5)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHannSlow(trSign*region.traces,num,'no'); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,6)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHPDer(trSign*region.traces,num,'no'); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

subplot(7,1,7)
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP_HF(trSign*region.traces,num,'no',trNew); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;

%there is a bug in calciumdxEvent_DetSingTrHannSlow.m that results in event offsets being out of bounds and truncated such that no. of onsets != no. of offsets and results in plotting error for hev2plottrace.m. Whatever this error is, it's been corrected in calciumdxEvent_FastAndSlow.m, which calls HannSlow but doesn't return in error (mismatched events have been checked). 