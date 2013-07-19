function data = batchFetchWaveProp(filelist,region)
%batchFetchWaveProperties Fetch calcium wave properties
%data = batchFetchWaveProp(filelist,region)
%
%Examples:
%data=batchFetchWaveProp({filename},region);
%data=batchFetchWaveProp(filelist);
%
%**USE**
%This function assumes that you have performed Calcium wave detection and the data is saved in all your region data structures that are being passed as filenames to this script.
%Options:
%Must provide one input:
%(1) table with desired filenames (space delimited txt file, with full filenames in first column)
%files.txt should have filenames in first column.
%can have extra columns with descriptor/factor information for the file. This will be the rowinfo that is attached to each measure observation in the following script.
%filelist = readtext('files.txt',' '); %grab readtext.m file script from matlab central
%or
%(2) a single filename (filename of your region .mat file) as a cell array, i.e.  {filename}
%filelist={filename}; %can pass just a single filename and a single already loaded region structure, if only getting values for a single file.
%Output:
%% When finished, convert 'data' table to string with the following code-- because matlab won't copy the contents of a mixed cell array correctly to other programs in the OS.
%tmp=data;
%for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
%tmp2=tmp';
%txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:})  %copy this output. **There will always be an extra column of tabs at end.
%%Change to the desired file name in the following line...
%dlmwrite('dWaveProps.txt',txt,'delimiter','','newline','unix');  %don't need this just copy output to console from above
%type('dWaveProps.txt');
%sprintf(['filename' '\t' 'ncells' '\t' 'pact' '\t' 'psynch' '\t' 'pgdp' '\t' 'psch' '\t' 'pother' '\t' 'psink' '\t' 'm_freq(Hz)' '\t' 'se_freq' '\t' 'm_dur(s)' '\t' 'se_dur' '\t' 'p.corr' '\t' 'n.pairs' '\t' 'pp.corr'...
%        '\r' filename '\t' num2str(ncells) '\t' num2str(pact) '\t' num2str(psynch) '\t' num2str(pgdp) '\t' num2str(psch) '\t' num2str(pother) '\t' num2str(psink) '\t' num2str(m_freq) '\t' num2str(se_freq) '\t' num2str(m_dur) '\t' num2str(se_dur) '\t' p_corr '\t' n_pairs '\t' p_pairs_corr])
%--------------------------------------------------------------------------
%Versions:
%James Ackman, 1/19/2011
%updated to work by location, 5/12/2011
%See also:
% batchFetchCalciumEventProps, calciumdxDetectWaves, calciumdxDetectWavesRefine, getWaveCentroids, getWaveSizeDistance, getActiveFraction, getWaveSpeeds, getWaveDirections

loadfile = 1;
if nargin > 1, loadfile = 0; end
% global region;
% results={};
results={'filename' 'region.name' 'wavenumber' 'nrois' 'roi.height.px' 'roi.width.px' 'waveonset.fr' 'waveoffset.fr' 'wavepeak.fr' 'nwaves' 'actvfraction' 'waveactvfraction' 'wavefreq.hz' 'waveisi.s' 'wavesize.um2' 'wavedist.um' 'wavespeed.umpersec' 'wavedir.degs'};

fnms=filelist(:,1);
for j=1:numel(fnms)
    
    if loadfile > 0
        matfile=load(fnms{j});
        region=matfile.region;
    end
    [pathstr, name, ext] = fileparts(fnms{j});
    rowinfo = {[name ext]};  %2011-07-11 jba
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
    
    
    locationMarkers = unique(region.location);
    regionsAll = splitRegion(region);
    for locationIndex = locationMarkers
        tmpregion = regionsAll{locationIndex}.region;
        try
            output=myWaveProps(tmpregion,rowinfo,locationIndex);
            results = [results; output;];
        catch exception
            %         if isempty(region.wavedata{locationIndex}.waveonsets)
            %             disp(['Region waveonsets is empty. Nothing to analyse for locationIndex = ' num2str(locationIndex)])
            %         else
            %             throw(exception);
            %         end
            rethrow(exception)
        end
    end
    h = waitbar(j/numel(fnms));
end
data=results;
assignin('base', 'data',data)
close(h)

tmp=data;
for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=tmp';
txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:})  %copy this output. **There will always be an extra column of tabs at end.
%Change to the desired file name in the following line...
dlmwrite('dWaveProps.txt',txt,'delimiter','','newline','unix');  %don't need this just copy output to console from above
end


function output = myWaveProps(region,rowinfo,locationIndex)
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
ROIsize = mean(strel_sz,1);

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
        if strncmp('2P',region.exptype,2)
        wavespeeds(i)=mean(consecutivespeeds);   %Default for 2P movies
        else
        wavespeeds(i)=median(consecutivespeeds); %if consecutive speeds not gaussian distributed, median can be more representative. Default for CCD movies
        end
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

%{
sprintf(['filename' '\t' 'nrois' '\t' 'nwaves' '\t' 'waveactvfraction' '\t' 'wavefreq.hz' '\t' 'waveisi.s' '\t' 'wavesize.um2' '\t' 'wavedist.um' '\t' 'wavespeed.umpersec' '\t' 'wavedir.degs'...
        '\r' filename '\t' num2str(nrois) '\t' num2str(nwaves) '\t' num2str(waveactvfraction) '\t' num2str(wavefreq) '\t' num2str(waveisi) '\t' num2str(wavesize) '\t' num2str(wavedist) '\t' num2str(wavespeed) '\t' num2str(wavedir)])
%}
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
