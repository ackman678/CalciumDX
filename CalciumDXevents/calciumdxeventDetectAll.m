%setup default detection params-------------------
if exist('prefs','var')
	parameterArray = setupDetectionPrefs2Params(prefs);	
else
	prefs = setupDetectionPreferences(region);
	parameterArray = setupDetectionPrefs2Params(prefsAll(1));
	detectionType = prefsAll(1).name;
end

NpopupDetect = get(popupDetect,'value');
detectorName = popupDetectList{NpopupDetect};
function_handle = str2func(detectorName);
hbar = waitbar(0,'Please wait...'); 
spk = zeros(size(nt));
dec = zeros(size(nt));
sz=size(region.traces);
for c = 1:sz(1) 
%     [s d] = calciumdxdettrialWaves(trSign*nt(c,:),region);
	%[s d] = feval(function_handle,trSign*nt(num,:),region,detectionType,'hannfilterorder',hannfilterorder,'sd',sd,'sd2',sd2,'sd3',sd3,'nonfilt',nonfilt,'hipass',hipass,'block_size',block_size,'start_baseline',start_baseline,'end_baseline',end_baseline,'maxOffsetTime',maxOffsetTime,'windowAverage',windowAverage,'baselineAverage',baselineAverage,'blockSizeMultiplier',blockSizeMultiplier,'slidingWinStartFrame',slidingWinStartFrame);
	[s d] = feval(function_handle,trSign*nt(num,:),region,detectionType,'parameterArray',parameterArray);

    if rem(c,10) == 0
        waitbar(c/sz(1),hbar); 
    end 
    spk(c,s) = 1; 
    dec(c,d) = 1; 
    if region.transients(1,c) == 1 && sum(spk(c,:)) > 0
        region.transients(1,c) = 4;
    end 
end
region.detectorname=['unsupervised calciumdxdettrial, ' detectionType ', ' parameterArray];
close(hbar)
hevPlotTrace
