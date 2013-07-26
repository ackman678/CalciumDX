function batchFetchCorrProps(filelist,region,datasetSelector)
%batchFetchCorrProps Batch Fetch Correlated pairs properties for each movie
%
%Examples:
%batchFetchCorrProps(filelist)
%batchFetchCorrProps({filename},region)
% 
%**USE**
%If you have run fetchCorrPairs.m  on your mat files and you have the 
%corr pairs data stored at region.userdata.corr_pairs,  then this script is
%very useful for calculating the no. of cells, percent cells correlated, 
%no. of pairs, and percent pairs correlated.  
%
%Must provide one input:
%(1) table with desired filenames (space delimited txt file, with full filenames in first column). Use 'readtext.m' from matlabcentral.
%filelist = readtext('files.txt',' ');
%Output:
%%Try: type('dCorrProps.txt');
%
%Versions:
%2013-07-18 15:40:37 James Ackman. Based on original batchFetchCorrPairs from 2007.09.28 by JBA.
%See also:
% fetchCorrPairs, batchFetchCorrPairs, batchFetchCalciumEventProps

if nargin < 3 || isempty(datasetSelector), datasetSelector = 1; end

if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end

matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
datafilename = fullfile(matlabUserPath,'dCorrProps.txt');
setupHeaders = exist(datafilename,'file');

if setupHeaders < 1
%setup headers for data set---------------------------------------------------------------
	if size(filelist,1) > 1 && size(filelist,2) > 1
	results2={'filename' 'matlab.filename' 'nCells' 'pCorr' 'nPairs' 'pPairsCorr'};
	else
	results2={'matlab.filename' 'nCells' 'pCorr' 'nPairs' 'pPairsCorr'};
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


    disp('--------------------------------------------------------------------')
    myCorrProps(region,rowinfo,datafilename);
    h = waitbar(j/numel(fnms));
end
close(h)






%---------------------------------
function myCorrProps(region,rowinfo,datafilename)

%locationMarkers = unique(region.location);  %TODO: for intersection of location indices with pairs indices
%all_s ={}; %jba

nCellsAll = size(region.traces,1);

%{
%--------------------------------------------------------------------------
%This section is for setting up all_s if you wanna do corr on the gdp, sch
%subsets of cells. Comment or uncomment if you wanna do these corrs.
if isfield(region,'userdata') %jba
	if isfield(region.userdata,'schmutzon') %jba
		%if isfield(region.userdata,'schmutzr') %added by jba
		all_s{2} = all_s{1}; % non-SCH cells
		all_s{2}(region.userdata.schmutzr,:) = [];
		all_s{3} = zeros(length(region.userdata.schmutzon),size(region.traces,2)); % SCH cells
		for i = 1:length(region.userdata.schmutzon)
			all_s{3}(i,region.userdata.schmutzon{i})=1;
		end
		str = {'all','non_sch','sch'};
	end %jba
else
	str = {'all'};
end %jba
%--------------------------------------------------------------------------
%}

output= {};   
%for ii = 1:length(all_s) %jba
    pairs = region.userdata.corr_pairs{1};
	%TODO: insert intersect location area indices with unique pairs indices here to give unique reshaped pair listing. May have to do combinations? Or just focus on inter regional values for now.
	nCells = nCellsAll;
    n_cells = nCells;
    p_corr = 100*length(unique(reshape(pairs,1,prod(size(pairs)))))/nCells;
    n_pairs = nCells*(nCells-1)/2;
    p_pairs_corr = 100*size(pairs,1)/n_pairs;
    output = [output; [rowinfo {n_cells p_corr n_pairs p_pairs_corr}];];
%end
fid = fopen(datafilename,'a');
for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=output';
fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);
