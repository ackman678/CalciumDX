%calciumdxeventDetect
NpopupDetect = get(popupDetect,'value');
detectorName = popupDetectList{NpopupDetect};
function_handle = str2func(detectorName);
spk(num,:) = 0;
dec(num,:) = 0;
% [s d] = calciumdxdettrialWaves(trSign*nt(num,:),region);
[s d] = feval(function_handle,trSign*nt(num,:),region);

spk(num,s) = 1;
dec(num,d) = 1;
if region.transients(1,num) == 1 && sum(spk(num,:)) > 0
    region.transients(1,num) = 4;
end
hevPlotTrace