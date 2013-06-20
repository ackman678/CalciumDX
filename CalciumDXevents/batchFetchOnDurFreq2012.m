%2012-07-19-- based on batchFetchOnDurFreq by JBA from 2007-10-28.  Updated for used for retinal wave recordings
%28.10.2007-- James Ackman
% Batch Fetch Onsets, durations, amplitudes, frequencies, and event intervals from
%If you have run find_calciumdxcorrpairs.m  on your mat files and you have the 
%saved .mat files output from calciumdx.
%This script is very useful for fetching a formatted list of these values from a table of filenames, 
% with optional extra data next to each filename (such as age, slice, layer, drug, etc).  This optional
%extra info will be concatenated to the results so that you are left with 5 nice tables of all your values 
%for easy input into data analysis programs (like R) for further analysis and plotting.

%Must provide one input:
%(1) table with desired filenames (space delimited txt file, with full filenames in first column)

%filelist = readtext('files.txt',' '); %grab readtext.m file script from matlab central

%% When finished, convert 'data' table to string-- matlab won't copy the contents of a mixed cell array correctly.
%Must convert to string, transpose, format with sprintf, then write to file with dlmwrite.

%{
tmp=data.freq;
for i=1:numel(tmp); tmp{i} = num2str(tmp{i}); end
tmp2=tmp';
txt=sprintf([repmat('%s\t',1,size(tmp2,1)),'\n'],tmp2{:});  %you can copy this output or use dlmwrite
%Change to the desired file name in the following line...
dlmwrite('dFreqs.txt',txt,'delimiter','','newline','unix');  %don't need this just copy output to console from above

type('data.txt');

%}


function data = batchFetchOnDurFreq2012(filelist)
%global region;


    
results.onsets={}; 
results.dur={}; 
results.freq={}; 
results.ampl={};
results.intvls={};

fnms=filelist(:,1);
    for j=1:numel(fnms)
    fnm=fnms{j};
    matfile=load(fnm);
    region=matfile.region;
    rowinfo = filelist(j,:);
    
    sprintf(fnm)
    
    ons=myOnsets(region,rowinfo);
    durs=myDur(region,rowinfo);
    freqs=myFreq(region,rowinfo);
    ampls=myAmpl(region,rowinfo);
    intvls=myInts(region,rowinfo);
    
    results.onsets=[results.onsets; ons;]; 
    results.dur=[results.dur; durs;]; 
    results.freq=[results.freq; freqs;]; 
    results.ampl=[results.ampl; ampls;];
    results.intvls=[results.intvls; intvls;];
    h = waitbar(j/numel(fnms));
    end

    data=results;

assignin('base', 'data',data)
close(h)
end



%--------------------------------------------------------------------
%fetch signal duration times for ea cell of a movie
function durs = myDur(region,rowinfo)
tmpres=region.timeres;
dur1 = [];
celltype = {};
cellnumb = {};
%mxdur = zeros(1,length(region.contours));
%mndur = zeros(1,length(region.contours));
for d = 1:length(region.contours)
    newdata = tmpres*(region.offsets{d}-region.onsets{d});
    dur1 = [dur1 newdata];
    
    %the following lines are used to fetch the active cell type string for adding to the dataframe
    tmpcell = {};
    tmpcellnumb = {};
    if ~isempty(region.onsets{d})
            for z = 1:length(newdata)
                tmpcell{z} = 'other';
                tmpcellnumb{z} = d;
            end
            celltype = [celltype tmpcell];
            cellnumb = [cellnumb tmpcellnumb];
    end
end

%the following few lines concatenates the rowinfo cell array with our numerical results
dur1=dur1';
dur2=num2cell(dur1);
celltype=celltype';
cellnumb=cellnumb';
tmp = {};
for m = 1:numel(dur1);
  tmp = [tmp; rowinfo;];
end
    
dur = [tmp cellnumb celltype dur2];
durs=dur;
end


%--------------------------------------------------------------------
%fetch signal onset duration times for ea cell of a movie
function ons = myOnsets(region,rowinfo)
tmpres=region.timeres;
ondur = {};
celltype = {};
cellnumb = {};
%values=[]; %TESTING
for c = 1:length(region.contours)
d=region.onsets{c};
e=region.offsets{c};
cellnumb{1} = num2str(c);    
	if ~isempty(region.onsets{c})
	tmpcell = 'other';
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
			celltype{1} = tmpcell;
			if ~isempty(seq2)
			idx=seq2(1);
			idx1=idx-1;
			value=idx1*tmpres;
%			values = [values; value]; %TESTING
			value=num2cell(value);
			ondur=[ondur; [rowinfo cellnumb celltype value];];
			else
			value = num2cell(NaN);
			ondur=[ondur; [rowinfo cellnumb celltype value];];
			end
		end
	end
end
ons=ondur;
end



%--------------------------------------------------------------------
%fetch frequencies (Hz) for ea cell of a movie
function freqs = myFreq(region,rowinfo)
tmpres=region.timeres;
actv = [];
s = rast2mat(region.onsets,size(region.traces,2));
actv = sum(s,2);
actv = actv/(size(region.traces,2)*tmpres);

idx(:,1)= 1:numel(actv); %for including all cells
%idx = find(actv > 0);
%actv = actv(find(actv>0)); %if you wanna get rid of zero frequencies

cellnumb = num2cell(idx);  %the cell numbers for our dataframe
celltype = {};
%the active cell type for our dataframe
for y = 1:numel(idx)
    if ~isempty(region.onsets{idx(y)})
                celltype{y} = 'other';
    else
        celltype{y} = []; %for all cells, including non-active ones
    end
end


%the following few lines concatenates the rowinfo cell array with our numerical results
actv1=num2cell(actv);
tmp = {};
for m = 1:numel(actv);
  tmp = [tmp; rowinfo;];
end
celltype = celltype';
actv = [tmp cellnumb celltype actv1];

freqs=actv;
end



%-------------------------------------------------------------------
%fetch signal amplitudes for ea cell of a movie
function ampls = myAmpl(region,rowinfo)
ampz={};
celltype = {};
cellnumb = {};

for c = 1:length(region.contours)
d=region.onsets{c};
e=region.offsets{c};

if ~isempty(region.onsets{c})
                tmpcell = 'other';
end
cellnumb{1} = num2str(c);

if ~isempty(region.onsets{c})
for k=1:length(d)

a=(region.traces(c,d(k):e(k)));
b=(a-a(1))/a(1); %normalize to baseline

%b=abs(b);
%amp = -1*max(b); %reverse since negative changes represents the fura2 amplitudes
amp = max(b); %for OGB or GCaMP
celltype{1} = tmpcell;
amp=num2cell(amp);
ampz=[ampz; [rowinfo cellnumb celltype amp];];
end

end
end
ampls = ampz;
figure; hist(cell2mat(ampz(:,4)),20); title('a-a(1)/a(1);')
end



%--------------------------------------------------------------------
function intvls = myInts(region,rowinfo)
tmpres=region.timeres;
intvls={};

for c = 1:length(region.contours)
d=region.onsets{c};
e=region.offsets{c};

% if any(d==0)
%     id=find(d==0);
%     d(id)=1;
% end

    if ~isempty(region.onsets{c})
                tmpcell = 'other';
    tmpcellnumb = c;
    
    d1 = [d size(region.traces,2)];
    e1 = [0 e-d];
    ints=diff([0 d1]) - e1;
    ints=ints*tmpres;
    
    celltype = {};
    cellnumb = {};
    for z = 1:numel(ints)
        celltype{z} = tmpcell;
        cellnumb{z} = num2str(tmpcellnumb);
    end

    ints=ints';
    celltype=celltype';
    cellnumb=cellnumb';

    %the following few lines concatenates the rowinfo cell array with our numerical results
    ints1=num2cell(ints);
    tmp = {};
    for m = 1:numel(ints);
      tmp = [tmp; rowinfo;];
    end

    ints2 = [tmp cellnumb celltype ints1];
    intvls=[intvls; ints2;];
    end
end
end

