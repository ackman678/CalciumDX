%setup default detection params-------------------
detectionType = 'normal';
hannfilterorder = 2;
sd=2;  %number of std dev to make threshold for detection
sd2=1;   %keep at 1 sd if there are large artifacts present in the same direction as the signal that are bigger than the signals of interest.
sd3=1;  %no. of stddev above average whole trace change.  2 is fine unless large baseline shifts in traces or big pos artifacts.
nonfilt = 1; %true(1) or false(0). Determines whether or not to use the orignal, non-low pass filtered data for signal onset refinement around line 215. Can be set to nonfilt = 1 (true) for already filtered data, or data with low acquisition rates (low noise during signal rise times).
hipass = 'true';
%----setup defaults for block size and F0 baseline window size-------------
%old def for inmed data (~10Hz in vitro fura2-AM imaging)
block_size = 3; %delta spacing for spikes and length of window in which to look for local minima in points.
% start_baseline = 50; %baseline in pts 
% end_baseline = 5; %baseline in pts	
%block_size = round(3/region.timeres); %%delta spacing for spikes and length of window in which to look for local minima in points.
start_baseline = round(22/region.timeres); %22sec baseline. Converted to no. of frames here.
end_baseline = round(4.5/region.timeres); %4.5 sec end baseline.
maxOffsetTime = round(5/region.timeres); %in sec.  15fr at 10Hz or 1.5sec is the one originally used with this code

NpopupDetect = get(popupDetect,'value');
detectorName = popupDetectList{NpopupDetect};
function_handle = str2func(detectorName);
hbar = waitbar(0,'Please wait...'); 
spk = zeros(size(nt));
dec = zeros(size(nt));
sz=size(region.traces);
for c = 1:sz(1) 
%     [s d] = calciumdxdettrialWaves(trSign*nt(c,:),region);
    [s,d,sd,sd2,sd3,nonfilt,hannfilterorder] = feval(function_handle,trSign*nt(num,:),region,detectionType,'hannfilterorder',hannfilterorder,'sd',sd,'sd2',sd2,'sd3',sd3,'nonfilt',nonfilt,'hipass',hipass,'block_size',block_size,'start_baseline', start_baseline,'end_baseline', end_baseline,'maxOffsetTime',maxOffsetTime);
    if rem(c,10) == 0
        waitbar(c/sz(1),hbar); 
    end 
    spk(c,s) = 1; 
    dec(c,d) = 1; 
    if region.transients(1,c) == 1 && sum(spk(c,:)) > 0
        region.transients(1,c) = 4;
    end 
end
region.detectorname=['unsupervised calciumdxdettrial, sd=' num2str(sd) ', sd2=' num2str(sd2) ', sd3=' num2str(sd3) ', nonfilt=' num2str(nonfilt) ', HannOrder=' num2str(hannfilterorder)];
close(hbar)
hevPlotTrace
