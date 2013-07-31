%calciumdxeventDetect
%setup default detection params-------------------
if exist('prefs','var')
	parameterArray = setupDetectionPrefs2Params(prefs);	
else
	prefsAll = setupDetectionPreferences(region);
	parameterArray = setupDetectionPrefs2Params(prefsAll(1));
	detectionType = prefsAll(1).name;
end

%detectionType = 'normal';
%hannfilterorder = 2;
%sd=3;  %number of std dev to make threshold for detection
%sd2=1;   %keep at 1 sd if there are large artifacts present in the same direction as the signal that are bigger than the signals of interest.
%sd3=2;  %no. of stddev above average whole trace change.  2 is fine unless large baseline shifts in traces or big pos artifacts.
%nonfilt = 0; %true(1) or false(0). Determines whether or not to use the orignal, non-low pass filtered data for signal onset refinement around line 215. Can be set to nonfilt = 1 (true) for already filtered data, or data with low acquisition rates (low noise during signal rise times). Unless high low pass filter order needs to be used for smoothing the traces then make sure to set this parameter to 0
%hipass = 'false';
%slidingWinStartFrame = 2;
%blockSizeMultiplier = 2;
%block_size = 3; %delta spacing for spikes and length of window in which to look for local minima in points.
%start_baseline = round(50/region.timeres); %22sec baseline. Converted to no. of frames here.
%end_baseline = round(10/region.timeres); %4.5 sec end baseline.
%maxOffsetTime = round(20/region.timeres); %in sec.  15fr at 10Hz or 1.5sec is the one originally used with this code
%windowAverage = 'mean'; %string: 'mean' | 'median'    	 
%baselineAverage = 'all'; %string: 'window' | 'all' 


%{
makePlots = 'true';
sz = get(0,'screensize');
figure('Name','dxeventDetect','NumberTitle','off','position',[1 sz(4) sz(3) sz(4)]);

str = {['hannfilterorder: ' num2str(hannfilterorder)], ['sd: ' num2str(sd)], ['sd2: ' num2str(sd2)], ['sd3: ' num2str(sd3)], ... 
['nonfilt: ' num2str(nonfilt)], ['hipass: ' num2str(hipass)], ['block_size: ' num2str(block_size)], ['start_baseline: ' num2str(start_baseline)], ...
['end_baseline: ' num2str(end_baseline)], ['maxoffsettime: ' num2str(maxOffsetTime)], ['windowAverage: ' windowAverage], ['baselineAverage: ' baselineAverage], ...
['slidingWinStartFrame: ' num2str(slidingWinStartFrame)], ['blockSizeMultiplier: ' num2str(blockSizeMultiplier)]};
%annotation('textbox', [.2 .4, .1, .1], 'String', str);
annotation('textbox', [.8 .25, .1, .1], 'String', str, 'Interpreter', 'none');

currNum = num;
allNum = [1 4 6 8 12 16 45 46 49 757 814];
for i = 1:length(allNum)
subplot(3,4,i) 
num=allNum(i);
calciumdxdettrial(trSign*nt(num,:),region,detectionType,'hannfilterorder',hannfilterorder,'sd',sd,'sd2',sd2,'sd3',sd3,'nonfilt',nonfilt,'hipass',hipass,'block_size',block_size,'start_baseline',start_baseline,'end_baseline',end_baseline,'maxOffsetTime',maxOffsetTime,'windowAverage',windowAverage,'baselineAverage',baselineAverage,'blockSizeMultiplier',blockSizeMultiplier,'slidingWinStartFrame',slidingWinStartFrame,'makePlots',makePlots);
title(['cell no. ' num2str(num)])
zoom xon
axis tight
end
%}




NpopupDetect = get(popupDetect,'value');
detectorName = popupDetectList{NpopupDetect};
function_handle = str2func(detectorName);
spk(num,:) = 0;
dec(num,:) = 0;
% [s d] = calciumdxdettrialWaves(trSign*nt(num,:),region);
%[s d] = feval(function_handle,trSign*nt(num,:),region,detectionType,'hannfilterorder',hannfilterorder,'sd',sd,'sd2',sd2,'sd3',sd3,'nonfilt',nonfilt,'hipass',hipass,'block_size',block_size,'start_baseline',start_baseline,'end_baseline',end_baseline,'maxOffsetTime',maxOffsetTime,'windowAverage',windowAverage,'baselineAverage',baselineAverage,'blockSizeMultiplier',blockSizeMultiplier,'slidingWinStartFrame',slidingWinStartFrame);
[s d] = feval(function_handle,trSign*nt(num,:),region,detectionType,'parameterArray',parameterArray);

spk(num,s) = 1;
dec(num,d) = 1;
if region.transients(1,num) == 1 && sum(spk(num,:)) > 0
    region.transients(1,num) = 4;
end
hevPlotTrace