function [results,responseArray,responseFreq,latencyArray] = myMakeMultiPETHplot(region,stimuliToPlot,bn,nstimuli,makePlots)
%myMakeMultiPETHplot - Fetch the firing latencies and frequencies for calcium imaging data in a time window around stimulus triggers
%James B. Ackman, July 25, 2011
%
%[results,responseArray,responseFreq,latencyArray] = myMakeMultiPETHplot(region,[],[],[],0)
%'region' is the calciumdxevents data structure to analyze
%'stimuliToPlot' is a numeric vector of which stimulus types to analyse. e.g. [1 2], where 1 is the left eye stimulus, and 2 is the right eye stimulus
%'bn' is a numeric vector of length 2 for the amount of time to analyse surrounding each stimulus in milliseconds. e.g. [-3000 5000] for -3 to 5 secs around a stim.
%'nstimuli' is a numeric vector representing which individual stimulus numbers from each stimuliToPlot to analyse. e.g. [1:3] to just analyse the 1st 3 stim (like if you want to avoid analysing the subsequent habituating responses). 
%'makePlots' is logical, defaulting to 1 for plotting figures


if nargin<2 || isempty(stimuliToPlot); stimuliToPlot=1:numel(region.stimuli); end
if nargin<3 || isempty(bn);      bn=[-300 800];      end
if nargin<4 || isempty(nstimuli); nstimuli=1:numel(region.stimuli{stimuliToPlot(1)}.stimulusParams); end
if nargin<5 || isempty(makePlots); makePlots = 1; end
if strcmp(nstimuli,'default'); 
nstimuliSwitch = 1;
else
nstimuliSwitch = 0;
end



locationMarkers = unique(region.location);
regionsAll = splitRegion(region);

locationIndex = 0;  %only needed since the below lines are commented out---
%the following commented out code is not needed, unless you don't want to look at all regions

for locationIndex = 1:numel(locationMarkers)
    tmpregion = regionsAll{locationMarkers(locationIndex)}.region;
    disp(region.name{locationMarkers(locationIndex)})
    try
        [Hist_all,xout,meanLatencies,nSpikes]=myPETH(tmpregion,stimuliToPlot,bn,nstimuli,region.name{locationMarkers(locationIndex)},makePlots);
        %         axeshandles{locationIndex}.ax = ax;
        results{locationIndex}.name = region.name{locationMarkers(locationIndex)};
        results{locationIndex}.Hist_all = Hist_all;
        results{locationIndex}.xout = xout;
        results{locationIndex}.meanLatencies = meanLatencies;
        results{locationIndex}.nSpikes = nSpikes;
%         text(0.1,0.1,region.name{locationMarkers(locationIndex)})
    catch exception
        rethrow(exception)
        %         if isempty(region.wavedata{locationIndex}.waveonsets)
        %             disp(['Region waveonsets is empty. Nothing to analyse for locationIndex = ' num2str(locationIndex)])
        %             startcounter(locationIndex) = startcounter(locationIndex) +1;
        %         else
        %             throw(exception);
        %         end
    end
end
%}

disp('sum all regions')
[Hist_all,xout,meanLatencies,nSpikes,responseArray,latencyArray]=myPETH(region,stimuliToPlot,bn,nstimuli,[region.name{locationMarkers}],makePlots);
results{locationIndex+1}.name = 'sum all regions';
results{locationIndex+1}.Hist_all = Hist_all;
results{locationIndex+1}.xout = xout;
results{locationIndex+1}.meanLatencies = meanLatencies;
results{locationIndex+1}.nSpikes = nSpikes;
for numStim = stimuliToPlot
	if nstimuliSwitch > 0; nstimuli=1:numel(region.stimuli{numStim}.stimulusParams); end
    responseFreq{numStim} = sum(responseArray{numStim},2);
    responseFreq{numStim} = responseFreq{numStim} ./ numel(nstimuli);
    if makePlots>0
    figure;
    handle_stim = [];
    colordata = repmat(responseFreq{numStim},1,3);
    colordata = mat2gray(colordata);
    [X, map] = gray2ind(colordata);
    %         map = jet(size(map,1));
    map = gray(size(map,1));
    map = flipud(map);
    colormap(map);
    %         colordata = 1-colordata;
    %         axis equal
    imagesize = size(region.image);
    xlim([0 imagesize(2)])
    ylim([0 imagesize(1)])
    set(gca,'ydir','reverse');
    %set(gca,'color',[0 0 0]);
    set(gca,'xtick',[],'ytick',[]);
    hold on
    for c = 1:size(region.contours,2)
        ind = X(c,1);
        ind = ind + 1;
        clr = map(ind,:);
        %            plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0 1 1]);
        handle_stim(c) = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),clr);
        %             set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
        set(handle_stim(c),'EdgeColor',[0.9 0.9 0.9]);  %outline
        %                         set(handle_stim(c),'FaceColor',mycolors(numStim,:));  %outline
    end
    colorbar('location','East')
    incr = imagesize(1)/6;
    title([region.stimuli{numStim}.description ', response fraction for ' num2str(numel(nstimuli)) ' stimuli'])
    text(0,incr,['max = ' num2str(max(responseFreq{numStim}))])
    text(0,2*incr,['min = ' num2str(min(responseFreq{numStim}))])
    text(0,3*incr,['mean = ' num2str(mean(responseFreq{numStim}))])
    end
end

%filename %roi no. %region.name %roi size %normalized xloca %normalized yloca %region.stimuli{numStim}.description %normalized responseFreq %absolutefiringFreq(dFreq) %meanLatency %meanAmpl %meanDur