function [region] = batchmakeStimParamsWaveonsets(filelist, region, appendStimuli)
%batchmakeStimParamsWaveonsets - Convert wave onsets and offsets to stimulus parameters
% Batch convert retinal wave onsets from a list of 'region' data structure files into stimulusParams and append to the list of any existing stimulusParams. Depends on makeStimParams to add and convert the frame indices to microsecond event markers. 
% By default the script appends the stimulus data to 'region.stimuli' and saves the region structure back to disk, so be aware. The relevant save line (~line 76) can be commented out for testing purposes with a single file. 
% TODO: make automatic stimlulus descriptor selection flexible and interactive between ~lines 50-70
%
% INPUTS: 
%	filelist - cell array of fullfile path strings.  Make a plaintxt file of local file names.  Use readtext.m from matlab central to read in filelist.
%		* e.g. %filelist = readtext('testfiles_2pt5x.txt',' '); %grab readtext.m file script from matlab central
%		* %filelist={filename}; %can also pass just a single filename and a single already loaded region structure, if only getting values for a single file.
%	region - struct, a valid region data structure, see calciumdx
%	appendStimuli - a string that must be either 'true' or 'false'.  Default is 'true' to append new stimuli, not write over existing ones. Only set to 'false' if you are sure you want to write over the existing stimuli data in your region data files. 
% 
% Examples: 
%	region=batchFetchStimResponseProps({filename}, region);
%	region=batchFetchStimResponseProps({fnm}, region, 'false');
%	batchmakeStimParamsWaveonsets(filelist, [], 'true');
%	batchmakeStimParamsWaveonsets(filelist);
%
% See also: makeStimParams, getStimParams, myFrameTriggerDetect, calciumdx, myBatchFilter
%
% Author: James B. Ackman 2/20/2012

if nargin< 3 || isempty(appendStimuli); appendStimuli = 'true'; end
if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end


fnms=filelist(:,1);
for j=1:numel(fnms)
    if loadfile > 0
        matfile=load(fnms{j});
        region=matfile.region;
    end
    
%     [pathstr, name, ext] = fileparts(fnms{j});
    %     rowinfo = {[name ext]};  %2011-07-11 jba
    %     rowinfo = filelist(j,:);
    sprintf(fnms{j})
    
    if strcmp(appendStimuli, 'false')
		region.stimuli = {};   %****------only if you want to delete old entries-----****
    end
    
    locationMarkers = unique(region.location);
    for locationIndices = locationMarkers
        if strcmp(region.name(locationIndices),'SC.R')
            desc1 = 'waveonsets.SC.R';
            desc2 = 'waveoffsets.SC.R';
        elseif strcmp(region.name(locationIndices),'SC.L')
            desc1 = 'waveonsets.SC.L';
            desc2 = 'waveoffsets.SC.L';
        elseif strcmp(region.name(locationIndices),'V1.R')
            desc1 = 'waveonsets.V1.R';
            desc2 = 'waveoffsets.V1.R';
        elseif strcmp(region.name(locationIndices),'V1.L')
            desc1 = 'waveonsets.V1.L';
            desc2 = 'waveoffsets.V1.L';
		elseif strcmp(region.name(locationIndices),'V2.R')
            desc1 = 'waveonsets.V2.R';
            desc2 = 'waveoffsets.V2.R';
        elseif strcmp(region.name(locationIndices),'V2.L')
            desc1 = 'waveonsets.V2.L';
            desc2 = 'waveoffsets.V2.L';
        else
            error('stimulus descriptor not recognized')
        end
        [region] = makeStimParams(region,region.wavedata{locationIndices}.waveonsets,desc1);
        [region] = makeStimParams(region,region.wavedata{locationIndices}.waveoffsets,desc2);
    end
    
    save(fnms{j},'region')
end