function region = calciumdxDetectWaves(region)
%calciumdxDetectWaves
%James Ackman, 1/6/2011
%need to have 'region' data structure loaded into workspace

%--wave detection script-----------------------------------------------------
%---Prepare input data in terms of a spike histogram------ This is done automatically below.
% data = mean(ntSpk,1);
% data(data > 0) = 1;

%---------------------------------------------
%--Build script to intersect nt trace image values with events-------------
%----must first do hipass filtering for baseline correction of region.traces to improve the intersected nt trace image of events (otherwise the colormap is scaled to the first waves, if there are strong baseline changes in mean fluorescence)---
Nyq=0.5*(1/region.timeres);
hipasscutoff = 0.005;
ntFilt=zeros(size(region.traces));

if 900 > size(region.traces,2)
    highfilterorder=round((size(region.traces,2)/3)) - 1;  %data must of length of greater than 3x filter order.
    
elseif size(region.traces,2) == 900
    highfilterorder=round((size(region.traces,2)/3)) - 2;
else
    highfilterorder = 300;
end

for i=1:size(region.traces,1)
    xf=filtfilt(fir1(highfilterorder,hipasscutoff,'high'),1,region.traces(i,:));
    ntFilt(i,:)=xf;
end
ntFilt=mat2gray(ntFilt);


%----START region.location based script------------------------------------

locationMarkers = unique(region.location);
% locationIndex = 1;
for locationIndex = locationMarkers
    f = find(region.location == locationIndex);
    numCell=numel(f);
%     spk = zeros(numel(f),size(region.traces,2));
%     dec = zeros(numel(f),size(region.traces,2));
    ntSpk = zeros(numCell,size(region.traces,2));
    matAttOnOff=zeros(numCell,size(region.traces,2));
for ind = 1:numCell
   if ~isempty(region.onsets{f(ind)})
      for d=1:length(region.onsets{f(ind)})
          ntSpk(ind,region.onsets{f(ind)}(d):region.offsets{f(ind)}(d)) = ntFilt(f(ind),region.onsets{f(ind)}(d):region.offsets{f(ind)}(d));
          matAttOnOff(ind,region.onsets{f(ind)}(d):region.offsets{f(ind)}(d)) = 1;
      end
   end
end

%{
%------Now intersect the new ntFilt values with event times for a thresholded event matrix coded by event amplitudes
% figure; imshow(mat2gray(ntFilt)); colormap(jet)
ntSpk=zeros(size(region.traces));
for c=1:size(region.traces,1)
    if ~isempty(region.onsets{c})
        for d=1:length(region.onsets{c})
            %            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = region.traces(c,region.onsets{c}(d):region.offsets{c}(d));
            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = ntFilt(c,region.onsets{c}(d):region.offsets{c}(d));
        end
    end
end
% figure; imshow(mat2gray(ntSpk)); colormap(flipud(gray))   %this flips the colormap so we have a white background
% figure; imshow(mat2gray(ntSpk)); colormap(jet)
%---------------------------------------------

    %use ntSpk creation method for onsets:offsets matrix.  delete the following matAttONoff lines... not neeeded
    

%------this gives the input spike histogram with different y values (fraction of cells). But same shape as the mean(ntSpk).
matAttOnOff=zeros(length(region.contours),length(region.traces));
matAttOn=zeros(length(region.contours),length(region.traces));
matAttOff=zeros(length(region.contours),length(region.traces));
numCell=length(region.contours);
for i=1:numCell
    ons=region.onsets{i};
    ofs=region.offsets{i};
    for j=1:length(ons)
        %         plot([ons(j) ofs(j)-1],[i i],'Linewidth',2)
        %         plot([ons(j)],[i ],'b.-')
        %         matAttOnOff(i,[ons(j):ofs(j)-1])=1;
        matAttOnOff(i,[ons(j):ofs(j)])=1;
        matAttOn(i,[ons(j)])=1;  %not uses.........
        matAttOff(i,[ofs(j)])=1;   %not used........
        %         covMat=xcov(region.traces(3,:)',sum(matAttOnOff)/numCell'
    end
end
%}
% data=sum(matAttOnOff)/numCell;
% % figure; plot(data)
% data(data > 0) = 1;
% % figure; plot(data)


%{
%------the following for loop is to smooth the square wave signal in case there are zero frames present in the middle of a wave peak due to incomplete transient detection. Should not be needed as long the auto detect parameters in calciumdxdettrial are not set to be overly stringent.
base_indices=find(data < 1);  %get zero frames
jitter_secs=10;
jitter=round(jitter_secs/region.timeres);
for i=base_indices
    if i > jitter && i <= numel(data)-jitter
        tmp=sum(data(i-jitter:i-1));
        tmp2=sum(data(i+1:i+jitter));
        if tmp > 0 && tmp2 > 0  %find any zero frames in the middle of a peak (which might result from incomplete calcium transient detection during a wave)
            data(i-jitter:i+jitter)=1;
        end
    end
end
% figure; plot(data) %in most instances, it will not have changed from above plot
%}


%-------the following will output the wave start times and end times and locations--------
% dx=diff(data);
% dx_up=dx;
% dx_up(dx_up < 0) = 0;
% % figure; plot(dx_up)
% dx_down=dx;
% dx_down(dx_down > 0) = 0;
% % figure; plot(-dx_down)

%waveonsets
% data=dx_up;
data=sum(matAttOnOff)/numCell;
data = myfilter(data,10);
deltaspacing=10; %10secs in between waves
mpd=round(deltaspacing/region.timeres); %wave onsets
[pks, idx] = findpeaks(data,'minpeakdistance',mpd,'threshold',0);
% idx1=idx1+1; %because the onsets are offset by one frame from diff
% figure;
% plot(data(1,:));
% hold on;
% plot(idx1,data(idx1),'ok');
% hold off;

figure;
a(1) = subplot(2,1,1);
plot(data(1,:));
hold on;
plot(idx,data(idx),'ok');
hold off;
a(2) = subplot(2,1,2);
imagesc(dfoverf(region.traces));
linkaxes(a,'x')

pks = [1 idx size(region.traces,2)];
minima = [];
for pkind = 1:(length(pks) - 1)
    datasegment = data(pks(pkind):pks(pkind+1));
    datasegmentMinima = find(datasegment == min(datasegment));
    datasegmentMinimum = fix(median(datasegmentMinima));
    minima = [minima (pks(pkind) + datasegmentMinimum-1)];
end
axes(a(1))
hold on
plot(minima,data(minima),'or');

wvonsets = [];
thresholdlevel = 0.05;
for pkind = 1:(length(minima) - 1)
    thresh = abs(data(minima(pkind)) - data(idx(pkind))) * thresholdlevel;
    datasegment = find(data(minima(pkind):idx(pkind)) > thresh+data(minima(pkind)));
    if isempty(datasegment)
        stn = 0;
    else
        stn = datasegment(1);
    end
    stn = stn + minima(pkind)-1;
    wvonsets = [wvonsets stn];
end
axes(a(1))
hold on
plot(wvonsets,data(wvonsets),'og');

wvoffsets = [];
thresholdlevel = 0.20;
for pkind = 1:(length(idx))
    thresh = abs(data(minima(pkind+1)) - data(idx(pkind))) * thresholdlevel;
    datasegment = find(data(idx(pkind):minima(pkind+1)) < thresh+data(minima(pkind+1)));
    if isempty(datasegment)
        stn = 0;
    else
        stn = datasegment(1);
    end
    stn = stn + idx(pkind)-1;
    wvoffsets = [wvoffsets stn];
end
axes(a(1))
hold on
plot(wvoffsets,data(wvoffsets),'ob');

idx1 = wvonsets;
idx2 = wvoffsets;

%--------------------------------------------------------------------------
%{
%waveoffsets
data=-dx_down;
mpd=round(deltaspacing/region.timeres); %wave offsets
[pks2, idx2] = findpeaks(data,'minpeakdistance',mpd);
% idx2=idx2+1; %because the offsets are offset by one frame from diff. Nope. They are fine from orig diff output.
figure;
plot(data(1,:));
hold on;
plot(idx2,data(idx2),'og');
hold off;
%}

%figure out if an offset was at last frame of movie (no. of onsets and offsets not equal)
if numel(idx1) ~= numel(idx2)
    button = questdlg('Number of onsets not equal to number of offsets (Event may be at end of movie). Set final offset to last frame of movie?');
    %     errordlg('Number of onsets not equal to number of offsets. Final offset set to last datapoint, in case wave was at end of movie.')
    if strcmp(button,'Yes')
        idx2=[idx2 size(region.traces,2)];
    end
end

%{
%wavepeaks
data = mean(ntSpk,1);
% deltaspacing=10; %10secs in between waves
% mpd=round(deltaspacing/region.timeres); %min inter wave interval converted to no. of frames
% [pks, idx] = findpeaks(data,'minpeakdistance',mpd);

pkidx=[];
pks = [];
for i=1:length(idx1)
    [C,I] = max(data(idx1(i):idx2(i)));
    mxidx = idx1(i) + I - 1;
    pkidx=[pkidx mxidx];
    pks = [pks C];
end
idx = pkidx;
%}
%--------------------------------------------------------------------------
region.wavedata{locationIndex}.waveonsets=idx1;
region.wavedata{locationIndex}.waveoffsets=idx2;
region.wavedata{locationIndex}.wavepeaks=idx;

%-------the following will output the wave peaks and locations in the histogram---
% figure;
% plot(data(1,:));
% hold on;
% plot(idx,data(idx),'or');
% plot(idx1,data(idx1),'ok');
% plot(idx2,data(idx2),'og');
% hold off;
end