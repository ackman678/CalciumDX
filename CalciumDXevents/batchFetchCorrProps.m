function batchFetchCorrProps(filelist,region,datasetSelector,locationMarkers)
%batchFetchCorrProps Batch Fetch Correlated pairs properties for each movie
%
%Examples:
%batchFetchCorrProps(filelist)  %batch over files
%batchFetchCorrProps({filename},region)  %single file
%batchFetchCorrProps({filename},region,[],[2 3]) %fetch props just for region.location markers 2 & 3 in addition to all
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
%or
%(2) a single filename (filename of your region .mat file) as a cell array, i.e.  {filename} with your region data structure loaded in workspace
%
%Options:
%region-- your region data structure loaded in workspace that you want to pass with your single {filename}
%datasetSelector-- a numeric integer n indicated which region.userdata.corr{n} dataset you want to fetch. Defaults to 1.
%locationMarkers-- a numeric vector of integers indicating which region.location indices/region.names you want include in the output. Defaults to all region.locations containing data (i.e. unique(region.location))
%
%Output:
%%Try: type('dCorrProps.txt');
%
%Versions:
%2013-07-18 15:40:37 James Ackman. Based on original batchFetchCorrPairs from 2007.09.28 by JBA.
%See also:
% fetchCorrPairs, batchFetchCorrPairs, batchFetchCalciumEventProps

if nargin < 3 || isempty(datasetSelector), datasetSelector = 1; end

if nargin < 4 || isempty(locationMarkers), 
%TODO: make secondary script for adding region.location and name and coords for additional cellType markers in addition to region type with interactive marking.
locationMarkers = unique(region.location);  %TODO: for intersection of location indices with pairs indices.  
end

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
	results2={'filename' 'matlab.filename' 'region.name' 'nCells' 'pCorr' 'nPairs' 'pPairsCorr'};
	else
	results2={'matlab.filename' 'region.name' 'nCells' 'pCorr' 'nPairs' 'pPairsCorr'};
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
    myCorrProps(region,rowinfo,datafilename,datasetSelector,locationMarkers);
    h = waitbar(j/numel(fnms));
end
close(h)



%---------------------------------
function myCorrProps(region,rowinfo,datafilename,datasetSelector,locationMarkers)
%Get corr properties for all cell pairs (default from fetchCorrPairs) as well as broken down by region.name and print table to txt file
%---Do all----------------
output= {}; 
str = {'all'};
nCellsAll = size(region.traces,1);  
pairs = region.userdata.corr{datasetSelector}.corr_pairs{1};
nCells = nCellsAll;
[p_corr,n_pairs,p_pairs_corr] = getCorrPairMetrics(nCells,pairs);
output = [output; [rowinfo str {nCells p_corr n_pairs p_pairs_corr}];];
printTable(str,nCells,p_corr,n_pairs,p_pairs_corr)

%---Do each locationMarker----------------
for locationIndex = locationMarkers
	str = {region.name{locationIndex}};
	cellIndices = find(region.location == locationIndex)';
	nCells = numel(cellIndices);
	%Now use ismember() to reshape the pairs listing based on the intersection with the cellIndices for the current locationMarker.
	%TODO: add combinations of all region.locations for printing of corr props
	Lia = ismember(pairs(:,1),cellIndices);  %intersect with 1st column of pairs indices
	pairs1 = pairs(Lia,:);
	Lia = ismember(pairs1(:,2),cellIndices); %intersect with 2nd column of pairs indices
	pairs1 = pairs1(Lia,:);
	[p_corr,n_pairs,p_pairs_corr] = getCorrPairMetrics(nCells,pairs1);
	output = [output; [rowinfo str {nCells p_corr n_pairs p_pairs_corr}];];
	printTable(str,nCells,p_corr,n_pairs,p_pairs_corr)
end

%---Now print output data to file----------
fid = fopen(datafilename,'a');
for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=output';
fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);



function [p_corr,n_pairs,p_pairs_corr] = getCorrPairMetrics(nCells,pairs)
p_corr = 100*length(unique(reshape(pairs,1,prod(size(pairs)))))/nCells;
n_pairs = nCells*(nCells-1)/2;
p_pairs_corr = 100*size(pairs,1)/n_pairs;



function printTable(str,nCells,p_corr,n_pairs,p_pairs_corr)
fprintf('\n');
fprintf(str{1});
fprintf('\n');
fprintf(['        Total number of cells: ' num2str(nCells) '\n']);
fprintf(['Percent cells in correlations: ' num2str(p_corr) '\n']);
fprintf(['        Total number of pairs: ' num2str(n_pairs) '\n']);
fprintf(['     Percent pairs correlated: ' num2str(p_pairs_corr) '\n']);




