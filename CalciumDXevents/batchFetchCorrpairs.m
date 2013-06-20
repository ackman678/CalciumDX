%28.09.2007-- James Ackman
% Batch Fetch Correlated pairs data
%If you have run find_calciumdxcorrpairs.m  on your mat files and you have the 
%corr pairs data stored at region.userdata.corr_pairs,  then this script is
%very useful for calculating the no. of cells, percent cells correlated, 
%no. of pairs, and percent pairs correlated.  

%The final output is a cell array titled
%'data' that will contain these values alongside factors for each set of 
%values (age, condition, cell correlation type (all cells, non-SCH cells, or SCH cells),etc)  

%Must provide one input:
%(1) table with desired filenames (space delimited txt file, with full filenames in first column)

%% When finished, convert 'data' table to string-- matlab won't copy the contents of a mixed cell array correctly.
%for i=1:numel(data); data{i} = num2str(data{i}); end


%data=data';
%txt=sprintf([repmat('%s\t',1,size(data,1)),'\n'],data{:})  %copy this output
%%dlmwrite('data.txt',txt,'');  %don't need this just copy output to console from above
%%type('data.txt');

%filelist = readtext('files.txt',' ');
function data = batchFetchCorrpairs(filelist)
global region;
results={};
    fnms=filelist(:,1);
    for j=1:numel(fnms)
    matfile=load(fnms{j});
    rowinfo = filelist(j,:);
    region=matfile.region;
    output=myCorr(region,rowinfo);
    results = [results; output;];
    h = waitbar(j/numel(fnms));
    end
data=results;
assignin('base', 'data',data)
close(h)
end
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
    p_corr = 100*length(unique(reshape(pairs,1,prod(size(pairs)))))/size(s,1);
    n_pairs = size(s,1)*(size(s,1)-1)/2;
    p_pairs_corr = 100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2);
    output = [output; [rowinfo str{ii} {n_cells p_corr n_pairs p_pairs_corr}];];
end
end
