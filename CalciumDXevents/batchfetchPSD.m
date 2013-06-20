%BatchFetchPSD
%JBA, Tuesday, October 21, 2008 3:55 PM
function [Hpsd, PxxAll,SampleFreqs] = batchfetchPSD(filelist,active)
%filelist is a list of calciumdx generated .mat files you want to process
%active is a switch (1 or 0) determining whether you want to limit the analysis to the active cell population
fnms=filelist(:,1);
PxxAll = [];
SampleFreqs = [];
figure();
    for j=1:numel(fnms)
    fnm=fnms{j};
    matfile=load(fnm);
    region=matfile.region;
    rowinfo = filelist(j,:);
    sprintf(fnm)
    
    nt = dfoverf(region.traces);
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
actvcells=find(sum(s,2)>0);
Fs = 1/region.timeres;
myCells = actvcells;
if (nargin > 1) && (active == 0)
myCells = setxor(actvcells, 1:size(nt,1)); %selects non-active cells
else
myCells = actvcells; %set of active cells
end
PxxM = myPSD(nt,myCells,Fs);
tmp = ones(1,numel(myCells));
tmp = tmp*Fs;
if size(PxxM,1) == 513
    PxxAll = [PxxAll PxxM];
    SampleFreqs = [SampleFreqs tmp];
else
    sprintf('%s',fnm,' --length of movie too long')
end
    
    
    
    h = waitbar(j/numel(fnms));
    end

    Fs = mean(SampleFreqs)
Hpsd = dspdata.psd(mean(PxxAll,2),'Fs',Fs);   %mean of the Power spectrum values
%Hpsd = dspdata.psd(mean(PxxAll,2),'Fs',8);
plot(Hpsd)  %plot the mean power spectrum
close(h)



%assignin('base', 'PxxAll',PxxAll)
%assignin('base', 'SampleFreqs',SampleFreqs)
end