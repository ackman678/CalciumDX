function [data, responseArray, latencyArray, results2] = batchFetchStimResponseProps(filelist,region,stimuliIndices,bn,nstimuliSelect,limitToWaves,limitToEventDataset,makePlots)
%batchFetchStimResponseProps - 	A wrapper and output generator for getting the firing latencies and frequencies for calcium imaging data in a time window around stimulus triggers
%James Ackman, 1/19/2011
%updated to work by location, JBA 5/12/2011
%updated to have to automatically save data in tab delimited output files, JBA 1/2012
%updated with event-based and cell-based dataset selection, 7/24/2012 JBA
%updated doc header JBA 5/1/2013
%Examples:
% >> data=batchFetchStimResponseProps({filename},region);
% >> data=batchFetchStimResponseProps(filelist);
% >> batchFetchStimResponseProps(filelist,[],[],[-5000 10000],[]);
% >> batchFetchStimResponseProps(filelist,[],'default',[-15000 30000],[]);
%
%**USE**
%Must provide one input:
%
%(1) table with desired filenames (space delimited txt file, with full filenames in first column)
%files.txt should have filenames in first column.
%can have extra columns with descriptor/factor information for the file. This will be the rowinfo that is attached to each measure observation in the following script.
%filelist = readtext('files.txt',' '); %grab readtext.m file script from matlab central
%or
%(2) a single filename (filename of your region .mat file) as a cell array, i.e.  {filename}
%
%Options:
%filelist={filename}; %can pass just a single filename and a single already loaded region structure, if only getting values for a single file.
%region is a 'region' data structure saved by calciumdx and calciumdxcalciumdexran
%stimuliIndices is either 'default' or a vector of stimulus type indices, i.e. [1 2].  This is which stimulus indices to analyse, i.e. if region.stimuli is a structure length of 2, then you have two types of stimuli and the stimulus indices are 1 and 2, or [1 2].   Default is 'default' which is for all stimulus types found in the region.stimulus structure.
%bn is a vector of length 2 for the time surrounding your stimuli in milliseconds you want to analyse. Default is [-300 800]
%nstimuliSelect is either 'default' or a numeric indicating of individual stimulus indices for each stimulus, i.e. [1 2 3 4 5], to analyze just the first 5 stimuli for each of stimulusIndices. Default is all stimuli for each type of stimulus.
%limitToWaves is 'true' or 'false'.  Default is 'false'.  This indicates if you want the analysed event frequencies and latencies in response to stimuli limited to those falling within wave periods, as defined by calciumdxDetectWaves.m
%limitToEventDataset is 'true' or 'false'.  Default is 'true'.  This indicates if you want to limit the output dataset to just an Event based dataset (bigger) vs a cell based dataset (smaller) with mean event latencies and props by row for ea cell.
%makePlots is 0 (false) or 1 (true) to indicate whether or not to make figure plots of mean triggered calcium event waveforms.
%
%Output:
%Right now this function will automatically write to two tab-delimited txt file outputs, a 'cell' based dataset 'dResponseFreq.txt' and an event based dataset 'dResponseFreqEvents.txt'.
%And these outputs will be appended if the file already exists. But there will be an extra copy of the column names in this case, so the file will have to be cleaned up in a text editor afterwards.
%
% See also myMakeMultiPETHplot, myPETH, getPETH, myFrameTriggerDetect, getStimParams, batchmakeStimParamsWaveonsets, makeStimParams, calciumdx, calciumdxevents, calciumdxDetectWaves

%% When finished, convert 'data' table to string with the following code-- because matlab won't copy the contents of a mixed cell array correctly to other programs in the OS.
%{
%Orig method of saving mixed cell array output, no longer necessary...
tmp=data;
for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=tmp';
txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:})  %copy this output. **There will always be an extra column of tabs at end.
%Change to the desired file name in the following line...
dlmwrite('dWaveProps.txt',txt,'delimiter','','newline','unix');  %don't need this just copy output to console from above

type('dResponseFreq.txt');
%}

%-----------------------------------------------------------------------------------------
%- Set up options and default parameters
%-----------------------------------------------------------------------------------------

if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end

if nargin<3 || isempty(stimuliIndices); stimuliIndices = 'default'; end
if nargin<4 || isempty(bn);  bn=[-300 800];      end
%if nargin<5 || isempty(nstimuli); nstimuli=1:numel(region.stimuli{stimuliToPlot(1)}.stimulusParams); end
if nargin<5 || isempty(nstimuliSelect); nstimuliSelect='default'; end
if nargin<6 || isempty(limitToWaves); limitToWaves = 'false'; end
if nargin<7 || isempty(limitToEventDataset); limitToEventDataset = 'true'; end
if nargin<8 || isempty(makePlots); makePlots = 0; end




% if nargin > 1, loadfile = 0; end
% global region;
% results={};
% results={'filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency' 'meanAmpl' 'meanDur' 'wavepeak.fr' 'nwaves' 'actvfraction' 'waveactvfraction' 'wavefreq.hz' 'waveisi.s' 'wavesize.um2' 'wavedist.um' 'wavespeed.umpersec' 'wavedir.degs'};

if size(filelist,1) > 1 && size(filelist,2) > 1
    results={'filename' 'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency'};
    results2={'filename' 'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'eventOnset.fr' 'meanLatency' 'eventNumber'};
else
    results={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency'};
    results2={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'eventOnset.fr' 'meanLatency' 'eventNumber'};
end
%filename %roi no. %region.name %roi size %normalized xloca %normalized yloca %region.stimuli{numStim}.description %normalized responseFreq %absolutefiringFreq(dFreq) %meanLatency %meanAmpl %meanDur

if strcmp(limitToEventDataset,'false')
    tmp=results;
    fid = fopen('dResponseFreq.txt','a');
    for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
    tmp2=tmp';
    fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
    fclose(fid);
end

tmp=results2;
fid = fopen('dResponseFreqEvents.txt','a');
for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=tmp';
fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);


fnms = filelist(:,1);
if size(filelist,1) > 1 && size(filelist,2) > 1
    fnms2 = filelist(:,2);
end
for j=1:numel(fnms)
    
    if loadfile > 0
        matfile=load(fnms{j});
        region=matfile.region;
    end
    
    
    if size(filelist,1) > 1 && size(filelist,2) > 1
        [pathstr, name, ext] = fileparts(fnms2{j});
        name1 = {[name ext]};  %2012-02-07 jba
        [pathstr, name, ext] = fileparts(fnms{j});
        name2 = {[name ext]};  %2012-02-07 jba
        rowinfo = [name1 name2];
    else
        [pathstr, name, ext] = fileparts(fnms{j});
        rowinfo = {[name ext]};  %2012-02-07 jba
    end
    
    %     rowinfo = filelist(j,:);
    sprintf(fnms{j})
    
    %    if nargin<3 || isempty(stimuliToPlot); stimuliToPlot=1:numel(region.stimuli); end
    if strcmp(stimuliIndices,'default');
        stimuliToPlot=1:numel(region.stimuli);
    else
        stimuliToPlot = stimuliIndices;
    end
    
    if strcmp(nstimuliSelect,'default');
        nstimuliSwitch = 1;
        nstimuli = 'default';
    else
        nstimuliSwitch = 0;
        nstimuli = nstimuliSelect;
    end
    
    
    %     %---Fix early 2011 2P JBA analysed files---
    %     region.wavedata{1}.waveonsets = region.waveonsets;
    %     region.wavedata{1}.waveoffsets = region.waveoffsets;
    %     region.wavedata{1}.wavepeaks = region.wavepeaks;
    %     region.wavedata{1}.wavecentroids = region.wavecentroids;
    %     region.wavedata{1}.waveprops = region.waveprops;
    %     region.wavedata{1}.waveprops.wavedirection.origin_coord = [0 0];
    %     region.wavedata{1}.name = region.brainarea;
    %     region.name = {}; region.name{1} = region.brainarea;
    %     %---END fix files--------------------------
    
    %---fix files in which onsets are located outside the detected waveframe intervals)--------
    if strcmp(limitToWaves,'true') && isfield(region,'wavedata')
        s=zeros(size(region.traces));
        for c=1:size(region.traces,1)
            locationIndex = region.location(c);
            if ~isempty(region.onsets{c})
                tempOnsets = [];
                tempOffsets = [];
                %             for d=1:length(region.onsets{c})
                for d = 1:numel(region.wavedata{locationIndex}.waveonsets);
                    %                 ind = find(region.wavedata{locationIndex}.waveonsets <= region.onsets{c}(d));
                    ind = find(region.onsets{c} >= region.wavedata{locationIndex}.waveonsets(d) & region.onsets{c}<= region.wavedata{locationIndex}.waveoffsets(d));
                    if ~isempty(ind)
                        %                     if region.onsets{c}(d) >= region.wavedata{locationIndex}.waveonsets(ind(end)) && region.onsets{c}(d) <= region.wavedata{locationIndex}.waveoffsets(ind(end))
                        s(c,region.onsets{c}(ind(1))) = 1;
                        %tempOnsets = [tempOnsets region.onsets{c}(ind(1))];
                        tempOnsets = [tempOnsets region.onsets{c}(ind)];
                        %tempOffsets = [tempOffsets region.offsets{c}(ind(1))];
                        tempOffsets = [tempOffsets region.offsets{c}(ind)];
                        %                     end
                    end
                end
                region.onsets{c} = tempOnsets;
                region.offsets{c} = tempOffsets;
            end
        end
    end
    %----END fix files------------------------------------------------------------------
    
    
    
    [histResults,responseArray,responseFreq,latencyArray] = myMakeMultiPETHplot(region,stimuliToPlot,bn,nstimuli,makePlots);
    
    %     locationMarkers = unique(region.location);
    %     regionsAll = splitRegion(region);
    %     for locationIndex = locationMarkers
    %         tmpregion = regionsAll{locationIndex}.region;
    %         try
    %             output=myWaveProps(tmpregion,rowinfo,locationIndex);
    %             results = [results; output;];
    %         catch exception
    %             %         if isempty(region.wavedata{locationIndex}.waveonsets)
    %             %             disp(['Region waveonsets is empty. Nothing to analyse for locationIndex = ' num2str(locationIndex)])
    %             %         else
    %             %             throw(exception);
    %             %         end
    %             rethrow(exception)
    %         end
    %     end
    disp('--------------------------------------------------------------------')
    for numStim = stimuliToPlot
        if nstimuliSwitch > 0; nstimuli=1:numel(region.stimuli{numStim}.stimulusParams); end
        if strcmp(limitToEventDataset,'false')
            output=myWaveProps(region,rowinfo,[],numStim,nstimuli,responseFreq,latencyArray);
            %results = [results; output;];
            
            fid = fopen('dResponseFreq.txt','a');
            for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
            tmp2=output';
            fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
            fclose(fid);
            clear output
        end
        
        %%        output2=myEventProps(region,rowinfo,[],numStim,nstimuli,responseFreq,latencyArray);
        myEventProps(region,rowinfo,[],numStim,nstimuli,responseFreq,latencyArray);
        %       %results2 = [results2; output2;];
        %}
    end
    
    h = waitbar(j/numel(fnms));
end
data=results;
% assignin('base', 'data',data)
close(h)

%tic
%tmp=data;
%for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
%tmp2=tmp';
%txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:});  %copy this output. **There will always be an extra column of tabs at end.
%%Change to the desired file name in the following line...
%dlmwrite('dResponseFreq.txt',txt,'delimiter','','newline','unix');
%toc

% ex1 = {'a' 1 12 123; 'ab' 4 5 6; 'abc' 7 8 9}
% tic
% tmp=ex1;
% fid = fopen('mydata2.txt','a');
% for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
% tmp2=tmp';
% fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
% fclose(fid);
% toc

%{
tic
tmp=results2;
for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=tmp';
txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:});  %copy this output. **There will always be an extra column of tabs at end.
%Change to the desired file name in the following line...
dlmwrite('dResponseFreqEvents.txt',txt,'delimiter','','newline','unix');
toc
%}
end


function output = myWaveProps(region,rowinfo,locationIndex,numStim,nstimuli,responseFreq,latencyArray)
output= {};

% filelist='/Users/ackman/Data/2photon/100913/TSeries-09132010-1150-040_defaultwithFindPeaks.mat';

%number of ROIs-- single value repeat for ea wave in data.frame------------
nrois=size(region.traces,1);

%get mean ROI size---------------------------------------------------------
%--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
strel_sz = zeros(length(region.contours),2);
for c = 1:length(region.contours)
    % c=1;
    % f={};
    [szX,szY] = size(region.image);  %assuming szY is the largest dimension
    % szZ = size(region.traces,2);
    if mod(max([szY szX]),min([szY szX])) == 0
        rXY=szY/szX;
        szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
    else
        rXY = 1;  %assuming other rectanqular images (like from CCD camera) don't have evenly divisible number of row lines as column lines (unlike for half scannin with a laser scanning microscope)
    end
    
    ps = round(region.contours{c});
    ps=[ps(:,1) rXY*ps(:,2)];
    if rXY > 1
        idx=find(ps(:,2) == min(ps(:,2)));
        ps(idx,2)=min(ps(:,2))-rXY;
    end
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    % f{c} = sub2ind([szX, szY],fy,fx);
    strel_sz(c,:) = [numel(unique(fx)) numel(unique(fx))];
end
ROIsize = mean(strel_sz,1) + 1;
%{
if isfield(region,'wavedata')
	if ~isempty(region.wavedata{locationIndex}.waveonsets)
		
		%number of waves-- single value repeat for ea wave in data.frame-----------
		nwaves=length(region.wavedata{locationIndex}.waveonsets);
		
		%actv fraction-- single value repeat for ea wave in data.frame-------------
		spk = zeros(size(region.traces));
		for c = 1:size(spk,1)
			spk(c,region.onsets{c}) = 1;
		end
		actv = sum(spk,2);
		actvfraction = length(find(actv>0))/size(spk,1);
		
		%waveactv fraction-- a vector the length(nwaves). Attached wave no. to ea measure value in order, so that last wave gets NaN.
		waveactvfraction = region.wavedata{locationIndex}.waveprops.waveactvfraction;
		
		%wave frequency-- single value repeat for ea wave in data.frame------------
		wavefreq= length(region.wavedata{locationIndex}.waveonsets)/(size(region.traces,2).*region.timeres);
		
		%wave isi-- a vector the length(nwaves). Attached wave no. to ea measure value in order, so that last wave gets NaN.
		tmpres=region.timeres;
		d=region.wavedata{locationIndex}.waveonsets;
		e=region.wavedata{locationIndex}.waveoffsets;
		
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
		% waveisi= mean(ints);
		waveisi=[ints NaN];  %since no. of measure isis will equal no. of waveframes-1, we add NaN to last wave. This way each ISI is attached to a wave (can analyse, size of wave vs next ISI for example) except for the last wave.
		wavesizemicrons2 = [];
		%wave size-----------------------------------------------------------------
		for i=1:numel(region.wavedata{locationIndex}.waveprops.waveareapixels)
			wavesizemicrons2(i) = (sum(region.wavedata{locationIndex}.waveprops.waveareapixels{i}) * region.spaceres^2);
		end
		% wavesize= mean(wavesizemicrons2(~isnan(wavesizemicrons2)));
		wavesize=wavesizemicrons2;
		
		%wave distance-------------------------------------------------------------
		% wavedist=mean(region.wavedata{locationIndex}.waveprops.wavedistance_px(~isnan(region.wavedata{locationIndex}.waveprops.wavedistance_px)));
		wavedist=region.wavedata{locationIndex}.waveprops.wavedistance_px .* region.spaceres;
		wavespeeds = [];
		%wave speed----------------------------------------------------------------
		for i=1:length(region.wavedata{locationIndex}.waveprops.wavespeeds_umpersec)
			consecutivespeeds=region.wavedata{locationIndex}.waveprops.wavespeeds_umpersec{i};
			%     wavespeeds(i)=mean(consecutivespeeds);   %Default for 2P movies
			wavespeeds(i)=median(consecutivespeeds); %if consecutive speeds not gaussian distributed, median can be more representative. Default for CCD movies
		end
		% wavespeed=mean(wavespeeds(~isnan(wavespeeds)));
		wavespeed=wavespeeds;
		
		%wave direction------------------------------------------------------------
		% wavedir = mean(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(~isnan(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians))) * (180/pi);
		centr = round(centroid(region.coords{locationIndex}));
		rowdifference = abs(region.wavedata{locationIndex}.waveprops.wavedirection.origin_coord(1) - region.orientation.value(1));
		if rowdifference < centr(2)
			theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians .* -1; %this flips the coordinate angles in one hemisphere so that the final angles are directly comparable between both hemipheres in the anterior (180deg), medial (270deg), lateral(90deg), and posterior(0deg) directions.
		else
			theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians;
		end
		% wavedir=(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians) .* (180/pi);
		wavedir = theta .* (180/pi);
		
		%loop through each wave, to get the following values-----------------------
		for wfr=1:length(region.wavedata{locationIndex}.waveonsets)
			output = [output; [rowinfo region.name{locationIndex} {wfr nrois ROIsize(1) ROIsize(2) region.wavedata{locationIndex}.waveonsets(wfr) region.wavedata{locationIndex}.waveoffsets(wfr) region.wavedata{locationIndex}.wavepeaks(wfr) nwaves actvfraction waveactvfraction(wfr) wavefreq waveisi(wfr) wavesize(wfr) wavedist(wfr) wavespeed(wfr) wavedir(wfr)}];];
		end
		
	else
		nwaves = 0;
		output = [output; [rowinfo region.name{locationIndex} {NaN nrois ROIsize(1) ROIsize(2) NaN NaN NaN nwaves 0 NaN 0 NaN NaN NaN NaN NaN}];];
	end
end
%}
latencies = nanmean(latencyArray{numStim},2);
locationMarkers = unique(region.location);
actvResponseFraction = cell(size(region.name));
for j = locationMarkers
    idx = find(region.location == j);
    actvResponseFraction{j} = numel(find(responseFreq{numStim}(idx) > 0)) / numel(responseFreq{numStim}(idx));
end

%loop through all rois, to get the following values-----------------------
for roiIND = 1:length(region.contours)
    nrois = numel(find(region.location == region.location(roiIND)));
    locapx = centroid(region.contours{roiIND});
    xlocapx = locapx(1);
    ylocapx = locapx(2);
    mx = max(region.coords{region.location(roiIND)});
    mn = min(region.coords{region.location(roiIND)});
    if region.orientation.value(1) < mn(2)
        xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
        ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
        ylocanorm = abs(1-ylocanorm);
    else
        xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
        ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
    end
    absFreq = numel(region.onsets{roiIND})/(size(region.traces,2).*region.timeres);
    output = [output; [rowinfo region.name{region.location(roiIND)} {roiIND nrois ROIsize(1) ROIsize(2) xlocapx ylocapx xlocanorm ylocanorm} region.stimuli{numStim}.description {numel(nstimuli) actvResponseFraction{region.location(roiIND)} responseFreq{numStim}(roiIND) absFreq latencies(roiIND)}];];
end
% results={'filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency' 'meanAmpl' 'meanDur'};

%{
sprintf(['filename' '\t' 'nrois' '\t' 'nwaves' '\t' 'waveactvfraction' '\t' 'wavefreq.hz' '\t' 'waveisi.s' '\t' 'wavesize.um2' '\t' 'wavedist.um' '\t' 'wavespeed.umpersec' '\t' 'wavedir.degs'...
        '\r' filename '\t' num2str(nrois) '\t' num2str(nwaves) '\t' num2str(waveactvfraction) '\t' num2str(wavefreq) '\t' num2str(waveisi) '\t' num2str(wavesize) '\t' num2str(wavedist) '\t' num2str(wavespeed) '\t' num2str(wavedir)])
%}
end


function output = myEventProps(region,rowinfo,locationIndex,numStim,nstimuli,responseFreq,latencyArray)
output= {};

% filelist='/Users/ackman/Data/2photon/100913/TSeries-09132010-1150-040_defaultwithFindPeaks.mat';

%number of ROIs-- single value repeat for ea wave in data.frame------------
%nrois=size(region.traces,1);
nCells=length(region.contours);

%get mean ROI size---------------------------------------------------------
%--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
strel_sz = zeros(length(region.contours),2);
for c = 1:nCells
    % c=1;
    % f={};
    [szX,szY] = size(region.image);  %assuming szY is the largest dimension
    % szZ = size(region.traces,2);
    if mod(max([szY szX]),min([szY szX])) == 0
        rXY=szY/szX;
        szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
    else
        rXY = 1;  %assuming other rectangular images (like from CCD camera) don't have evenly divisible number of row lines as column lines (unlike for half scannin with a laser scanning microscope)
    end
    
    ps = round(region.contours{c});
    ps=[ps(:,1) rXY*ps(:,2)];
    if rXY > 1
        idx=find(ps(:,2) == min(ps(:,2)));
        ps(idx,2)=min(ps(:,2))-rXY;
    end
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    % f{c} = sub2ind([szX, szY],fy,fx);
    strel_sz(c,:) = [numel(unique(fx)) numel(unique(fx))];
end
ROIsize = mean(strel_sz,1) + 1;

%latencies = nanmean(latencyArray{numStim},2);
locationMarkers = unique(region.location);
actvResponseFraction = cell(size(region.name));
for j = locationMarkers
    idx = find(region.location == j);
    actvResponseFraction{j} = numel(find(responseFreq{numStim}(idx) > 0)) / numel(responseFreq{numStim}(idx));
end

%nCells=length(region.contours);

%loop through all rois, to get the following values-----------------------
for roiIND = 1:nCells
    nrois = numel(find(region.location == region.location(roiIND)));
    locapx = centroid(region.contours{roiIND});
    xlocapx = locapx(1);
    ylocapx = locapx(2);
    mx = max(region.coords{region.location(roiIND)});
    mn = min(region.coords{region.location(roiIND)});
    if region.orientation.value(1) < mn(2)
        xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
        ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
        ylocanorm = abs(1-ylocanorm);
    else
        xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
        ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
    end
    absFreq = numel(region.onsets{roiIND})/(size(region.traces,2).*region.timeres);
    
    nEvents = size(latencyArray{numStim},2);
    latencies = num2cell(latencyArray{numStim}(roiIND,:))';
    eventIndices = num2cell((1:nEvents))';
    %    disp(size(latencies))
    %    disp(size(eventIndices))
    
    tmpinfo = {};
    for eventIdx = 1:nEvents
        tmp = [rowinfo region.name{region.location(roiIND)} {roiIND nrois ROIsize(1) ROIsize(2) xlocapx ylocapx xlocanorm ylocanorm} region.stimuli{numStim}.description {numel(nstimuli) actvResponseFraction{region.location(roiIND)} responseFreq{numStim}(roiIND) absFreq region.stimuli{numStim}.stimulusParams{eventIdx}.frame_indices(1,1)}];
        tmpinfo = [tmpinfo; tmp;];
    end
    
    %    for eventIdx = 1:nEvents
    %		output = [output; [rowinfo region.name{region.location(roiIND)} {roiIND nrois ROIsize(1) ROIsize(2) xlocapx ylocapx xlocanorm ylocanorm} region.stimuli{numStim}.description {numel(nstimuli) actvResponseFraction{region.location(roiIND)} responseFreq{numStim}(roiIND) absFreq latencyArray{numStim}(roiIND,eventIdx) eventIdx}];];
    output = [tmpinfo latencies eventIndices];
    fid = fopen('dResponseFreqEvents.txt','a');
    for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
    tmp2=output';
    fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
    fclose(fid);
    %    end
    %    disp(roiIND)
    h2 = waitbar(roiIND/nCells);
end
close(h2)
end


%{
%---------------------------------
%fetch signal duration times for ea cell of a movie
function output = myCorr(region,rowinfo)
output= {};
all_s ={}; %jba
all_s{1} = zeros(size(region.traces)); % all cells
for i = 1:size(region.traces,1)
    all_s{1}(i,region.onsets{i})=1;
end

if isfield(region,'userdata') %jba
if isfield(region.userdata,'schmutzon') %jba
%if isfield(region.userdata,'schmutzr') %added by jba
all_s{2} = all_s{1}; % non-SCH cells
all_s{2}(region.userdata.schmutzr,:) = [];
if isfield(region.userdata,'schmutzon') %added by jba
all_s{3} = zeros(length(region.userdata.schmutzon),size(region.traces,2)); % SCH cells
for i = 1:length(region.userdata.schmutzon)
    all_s{3}(i,region.userdata.schmutzon{i})=1;
end
end %jba
end %jba
end %jba
str = {'all','non_sch','sch'};

for ii = 1:length(all_s) %jba
    s = all_s{ii};
    pairs = region.userdata.corr_pairs{ii};

    n_cells = size(s,1);
    p_corr = 100*length(unique(reshape(pairs,1,numel(pairs))))/size(s,1);
    n_pairs = size(s,1)*(size(s,1)-1)/2;
    p_pairs_corr = 100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2);
    output = [output; [rowinfo str{ii} {n_cells p_corr n_pairs p_pairs_corr}];];
end
end
%}