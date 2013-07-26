function batchFetchCorrPairs(filelist,region,datasetSelector)
%batchFetchCorrPairs Batch Fetch Correlated pairs data
%
%Examples:
%batchFetchCorrPairs(filelist)  %batch over files
%batchFetchCorrPairs({filename},region) %single file
%batchFetchCorrPairs({filename},region,2) %use 2nd corr dataset (region.userdata.corr{2})
% 
%**USE**
%If you have run fetchCorrPairs.m  on your mat files and you have the 
%corr pairs data stored at region.userdata.corr_pairs,  then this script is
%very useful for outputting a data file with a pair list with pvalues and inter-cell centroid distances. Right now it provides an undirected pair list, useful for constructing symmetric graphs.
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
%
%Output:
%%Try: type('dCorrPairs.txt');
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
    myCorr(region,rowinfo,datafilename,datasetSelector);
    h = waitbar(j/numel(fnms));
end
close(h)




%---------------------------------
function myCorr(region,rowinfo,datafilename,datasetSelector)
%Get corr data (cell list, pvalues, centroid-centroid physical distance) for all cell pairs and print table to txt file
output= {}; 
pairs = region.userdata.corr{datasetSelector}.corr_pairs{1};

%--setup roi height width ratio--------------
%The following is important for getting the distances right if the data pixel dimensions are not equivalent
%And below the scripts will assume wherever 'rXY' is used, that it is szX (m dimension) which must be scaled up.
%the following assumes that the modulus of raster scanned data is 0 (equally divisible image size) and that for CCD images the ratio of image dimensions is either equivalent or not equally divisible
[szX,szY] = size(region.image);  %assuming szY is the largest dimension and szX may or may not need to be scaled.
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;
end
%-- end setup roi height width ratio---------


for i = 1:size(pairs,1)
	cellA = pairs(i,1);
	cellB = pairs(i,2);
	pvalue = region.userdata.corr{datasetSelector}.pvalCorrMatrix(cellA,cellB);
	dist_px = getCellCellDistance(region,cellA,cellB,rXY);
	output = [output; [rowinfo {cellA cellB pvalue dist_px}];];
end

%---Now print output data to file----------  21.918719 seconds outside of for loop.
fid = fopen(datafilename,'a');
for i=1:numel(output); output{i} = num2str(output{i}); end  %this will be to 4 decimal points (defaut for 'format short'). Can switch to 'format long' before running this loop if need more precision.
tmp2=output';
fprintf(fid,[repmat('%s\t',1,size(tmp2,1)-1),'%s\n'],tmp2{:});
fclose(fid);


function dist_px = getCellCellDistance(region,cellA,cellB,rXY)
centr1 = centroid(region.contours{cellA});
centr2 = centroid(region.contours{cellB});
dist_px = sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2);	
