function batchFetchCorrPairs(filelist,region,datasetSelector)
%batchFetchCorrpairs Batch Fetch Correlated pairs data
%
%Examples:
%batchFetchCorrPairs(filelist)
%batchFetchCorrPairs({filename},region)
% 
%**USE**
%If you have run fetchCorrPairs.m  on your mat files and you have the 
%corr pairs data stored at region.userdata.corr_pairs,  then this script is
%very useful for outputting a data file with a pair list with pvalues and inter-cell centroid distances. Right now it provides an undirected pair list, useful for constructing symmetric graphs.
%Must provide one input:
%(1) table with desired filenames (space delimited txt file, with full filenames in first column)
%
%Output:
%% When finished, convert 'data' table to string-- matlab won't copy the contents of a mixed cell array correctly.
%for i=1:numel(data); data{i} = num2str(data{i}); end
%data=data';
%txt=sprintf([repmat('%s\t',1,size(data,1)),'\n'],data{:})  %copy this output
%%dlmwrite('data.txt',txt,'');  %don't need this just copy output to console from above
%%type('data.txt');
%filelist = readtext('files.txt',' ');
%
%Versions:
%2013-07-18 15:40:37 James Ackman. Based on original batchFetchCorrPairs from 2007.09.28 by JBA.
%See also:
% fetchCorrPairs, batchFetchCorrProps, batchFetchCalciumEventProps

if nargin < 3 || isempty(datasetSelector), datasetSelector = 1; end

if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end

matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
datafilename = fullfile(matlabUserPath,'dCorrPairs.txt');
setupHeaders = exist(datafilename,'file');

if setupHeaders < 1
%setup headers for data set---------------------------------------------------------------
	if size(filelist,1) > 1 && size(filelist,2) > 1
	results2={'filename' 'matlab.filename' 'cellA' 'cellB' 'pvalue' 'dist.px'};
	else
	results2={'matlab.filename' 'cellA' 'cellB' 'pvalue' 'dist.px' };
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
    myCorr(region,rowinfo,datafilename);
    h = waitbar(j/numel(fnms));
end
close(h)




%---------------------------------
function myCorr(region,rowinfo,datafilename)

nCellsAll = size(region.traces,1);

output= {};   
    pairs = region.userdata.corr_pairs{1};
	nCells = nCellsAll;
    n_cells = nCells;
    p_corr = 100*length(unique(reshape(pairs,1,prod(size(pairs)))))/nCells;
    n_pairs = nCells*(nCells-1)/2;
    p_pairs_corr = 100*size(pairs,1)/n_pairs;
    output = [output; [rowinfo {n_cells p_corr n_pairs p_pairs_corr}];];
repat tmpinfo pvalue distance cat tp pairs for output
fid = fopen(datafilename,'a');
for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=output';
fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);
