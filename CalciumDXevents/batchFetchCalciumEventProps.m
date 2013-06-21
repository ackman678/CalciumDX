function batchFetchCalciumEventProps(filelist,region,limitToWaves)
%batchFetchCalciumEventProps - 	A wrapper and output generator for getting the amplitudes, firing frequencies, event rise times and durations, normalized ROI spatial locations for calcium imaging data
%James Ackman, 2012-07-23
%2012-07-23-- based on batchFetchOnDurFreq and batchFetchOnDurFreq2012 by JBA from 2007-10-28 - 2012-07-23.  Updated for used for retinal wave recordings
%Examples:
% >> batchFetchCalciumEventProps(filelist);
% >> batchFetchCalciumEventProps({filename},region);
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
%stimuliIndices is either 'default' or a vector of stimulus type indices, i.e. [1 2].  This is which stimulus indices to analyse, i.e. if region.stimuli is a structure length of 2, then you have two types of stimuli and the stimulus indices are 1 and 2, or [1 2].   Default is 'default' which is for all stimulus types found in the region.stimulus structure.
%bn is a vector of length 2 for the time surrounding your stimuli in milliseconds you want to analyse. Default is [-300 800]
%nstimuliSelect is either 'default' or a numeric indicating of individual stimulus indices for each stimulus, i.e. [1 2 3 4 5], to analyze just the first 5 stimuli for each of stimulusIndices. Default is all stimuli for each type of stimulus.
%limitToWaves is 'true' or 'false'.  Default is 'false'.  This indicates if you want the analysed event frequencies and latencies in response to stimuli limited to those falling within wave periods, as defined by calciumdxDetectWaves.m
%
%Output:
%Right now this function will automatically write to a tab-delimited txt file outputs, a 'cell' based dataset 'dResponseFreq.txt' and an event based dataset 'dEvents.txt'.
%And these outputs will be appended if the file already exists. But there will be an extra copy of the column names in this case, so the file will have to be cleaned up in a text editor afterwards.
%
% See also batchFetchStimResponseProps, myMakeMultiPETHplot, myPETH, getPETH, myFrameTriggerDetect, getStimParams, batchmakeStimParamsWaveonsets, calciumdx, calciumdxevents, calciumdxDetectWaves

%Versions
%2013-06-21 14:34:09 updated new default dEventProps.txt file location to default matlab userpath folder

%-----------------------------------------------------------------------------------------
%- Set up options and default parameters
%-----------------------------------------------------------------------------------------

if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end

if nargin<3 || isempty(limitToWaves); limitToWaves = 'false'; end

matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
datafilename = fullfile(matlabUserPath,'dEventProps.txt');
setupHeaders = exist(datafilename,'file');

if setupHeaders < 1
%setup headers for data set---------------------------------------------------------------
	if size(filelist,1) > 1 && size(filelist,2) > 1
	results2={'filename' 'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm' 'freq.hz' 'intvls.s' 'onsets.s' 'durs.s' 'ampl.df'};
	else
	results2={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm' 'freq.hz' 'intvls.s' 'onsets.s' 'durs.s' 'ampl.df'};
	end
	%filename %roi no. %region.name %roi size %normalized xloca %normalized yloca %region.stimuli{numStim}.description %normalized responseFreq %absolutefiringFreq(dFreq) %meanLatency %meanAmpl %meanDur


	%open file for writing--------------------------------------------------------------------
	tmp=results2;
	fid = fopen(datafilename,'a');
	for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
	tmp2=tmp';
	fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
	fclose(fid);
end

%start loop through files-----------------------------------------------------------------
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

	%---fix files in which onsets are located outside the detected waveframe intervals)--------   this should be put inside openingFnc/wrapperFnc 
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

    disp('--------------------------------------------------------------------')
	myEventProps(region,rowinfo,datafilename);
    h = waitbar(j/numel(fnms));
end
%data=results;
close(h)
end



%-----------------------------------------------------------------------------------------
function output = myEventProps(region,rowinfo,datafilename)
output= {};
%number of ROIs-- single value repeat for ea wave in data.frame------------
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

fid = fopen(datafilename,'a');

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
    
    nEvents = numel(region.onsets{roiIND});
    if ~isempty(region.onsets{roiIND})
		if length(region.onsets{roiIND}) == length(region.offsets{roiIND})
		tmpinfo = {};  %output instead?
		for eventIdx = 1:nEvents  %loop through for nEvents and create data frame.  Ons, durs and ampls should all be same size and will be horz cat to this cell array
			tmp = [rowinfo region.name{region.location(roiIND)} {roiIND nrois ROIsize(1) ROIsize(2) xlocapx ylocapx xlocanorm ylocanorm absFreq}];
			tmpinfo = [tmpinfo; tmp;];  %output instead?	
		end    
		
		ons = myOnsets(region,roiIND);
		durs = myDur(region,roiIND);
		ampls = myAmpl(region,roiIND);
		intvls = myInts(region,roiIND);
	
		%fastest to save data on each iteration, rather than keep a growing cell array ('output') and save at end, cause the growing cell array gets very slow.
		output = [tmpinfo intvls ons durs ampls];    
		for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
		tmp2=output';
		fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
		
	%	output = [output; [tmpinfo intvls ons durs ampls]];
		end
		
	else
		tmpinfo = {};  %output instead?
		tmp = [rowinfo region.name{region.location(roiIND)} {roiIND nrois ROIsize(1) ROIsize(2) xlocapx ylocapx xlocanorm ylocanorm absFreq}];
		tmpinfo = [tmpinfo; tmp;];  %output instead?	
		ons = {NaN};
		durs = {NaN};
		ampls = {NaN};
		intvls = {NaN};
		%fastest to save data on each iteration, rather than keep a growing cell array ('output') and save at end, cause the growing cell array gets very slow.
		output = [tmpinfo intvls ons durs ampls];    
		for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
		tmp2=output';
		fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});			
	end
	h2 = waitbar(roiIND/nCells);
end
%fid = fopen(datafilename,'a');
%for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
%tmp2=output';
%fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);
close(h2)
end



%-----------------------------------------------------------------------------------------
function durs = myDur(region,roiIND)
%fetch signal duration times for ea cell of a movie
tmpres=region.timeres;
dur1 = tmpres*(region.offsets{roiIND}-region.onsets{roiIND});
durs=num2cell(dur1');
end



%-----------------------------------------------------------------------------------------
function ons = myOnsets(region,roiIND)
%fetch signal onset duration times for ea cell of a movie
tmpres=region.timeres;
values=[];
c = roiIND;
d=region.onsets{c};
e=region.offsets{c};
	if ~isempty(region.onsets{c})
		for k=1:length(d)
		%min1=min(region.traces(c,d(k):e(k)));
		%idx1=find(region.traces(c,d(k):e(k))==min1);
		a=(region.traces(c,d(k):e(k)));
		b=(a-a(1))/a(1);
		
%		b=abs(b);
		seq=b(1:find(b==max(b)));
		seq2=find(seq>=(max(b)/2) & seq~=0); %half amplitude rise time
		%seq=b(1:find(b==min(b)));
		%seq2=find(seq<=(min(b)/2) & seq~=0); %half amplitude rise time
			if ~isempty(seq2)
			idx=seq2(1);
			idx1=idx-1;
			value=idx1*tmpres;
			else
			value = NaN;
			end
			values = [values; value];
		end
	end
ons=num2cell(values);
end



%-----------------------------------------------------------------------------------------
function ampls = myAmpl(region,roiIND)
%fetch signal amplitudes for ea cell of a movie
ampz=[];
c = roiIND;
d=region.onsets{c};
e=region.offsets{c};
if ~isempty(region.onsets{c})
	for k=1:length(d)
	a=(region.traces(c,d(k):e(k)));
	b=(a-a(1))/a(1); %normalize to baseline
	%b=abs(b);
	%amp = -1*max(b); %reverse since negative changes represents the fura2 amplitudes
	amp = max(b); %for OGB or GCaMP
	ampz = [ampz; amp];
	end
end
ampls = num2cell(ampz);
%figure; hist(cell2mat(ampz(:,4)),20); title('a-a(1)/a(1);')
end



%-----------------------------------------------------------------------------------------
function intvls = myInts(region,roiIND)
%Interevent intervals
tmpres=region.timeres;
c = roiIND;
d=region.onsets{c};
e=region.offsets{c};
d1 = [d size(region.traces,2)];
e1 = [0 e-d];
ints=diff([0 d1]) - e1;
if d(1) ~= 1
	ints=ints(2:end);
end
if d(end) ~= size(region.traces,2)
	ints=ints(1:end-1);
end
ints=[ints*tmpres NaN];  %since no. of measured intervals will be nEvents-1, cat an 'NaN' at the end of the cell array
intvls=num2cell(ints');
end