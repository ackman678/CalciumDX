function [hf1, trHist]=myPlotRasterHistManualPeaks(fnm,region,stimulusParams,rmBaseline,dur,spk,dec)
%plot spike rasterplot with overlying histogram
%if dur = 'true', then a onset-duration rasterplot and histogram is generated
%if rmBaseline = 'true', then baseline drift correction is performed
%usage: myPlotRasterHistManualPeaks(fnm,region,[],[],'true')
%James Ackman, 1/6/2011
if nargin < 5, dur = 'false'; end
if nargin < 4, rmBaseline = 'false'; end
if nargin < 3, stimulusParams = []; end

%Network activity plots-----------------------------------------------------
%Left histogram plot of ROI activity-------------------------------------
% nt = filtfilt(fir1(10,0.1,'low'),1,region.traces);
% nt = dfoverf(region.traces);
nt = zeros(size(region.traces));
for c = 1:size(region.traces,1)
    nt(c,:) = dfoverf(region.traces(c,:))*100;
end


hf1(1)=figure();
% scrsize = get(0,'screensize');
% set(gcf,'position',[scrsize(3)/2-8.5/11*0.86*scrsize(4)/2 0.07*scrsize(4) 8.5/11*0.86*scrsize(4)*0.5 0.86*scrsize(4)]);
% set(gcf,'color',[1 1 1]);
% set(gcf,'PaperType','usletter');
% set(gcf,'PaperPositionMode','auto');

ax(1)=subplot(3,1,1);
%{
if strcmp(dur,'true')
    s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
else
    s = rast2mat(region.onsets,size(region.traces,2));
end
%}

%----
% spk=region.onsets;
% sz=size(region.traces,2);
% s = zeros(size(spk,2),sz);
% for c = 1:size(spk,2)
%    s(c,spk{c}) = 1;
% end
% %----

%----alternate hist method, not involving patch---------------

onsets = cell(1,size(spk,1));
offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    onsets{c} = find(spk(c,:)==1);
    offsets{c} = find(dec(c,:)==1);
end

region.onsets=onsets;
region.offsets=offsets;

matAttOnOff=zeros(length(region.onsets),length(region.traces));
matAttOn=zeros(length(region.onsets),length(region.traces));
matAttOff=zeros(length(region.onsets),length(region.traces));
numCell=length(region.onsets);
for i=1:numCell
        ons=region.onsets{i};
        ofs=region.offsets{i};
    for j=1:length(ons)
%         plot([ons(j) ofs(j)-1],[i i],'Linewidth',2)
%         plot([ons(j)],[i ],'b.-')
        matAttOnOff(i,[ons(j):ofs(j)-1])=1;
        matAttOn(i,[ons(j)])=1;
        matAttOff(i,[ofs(j)])=1;
%         covMat=xcov(region.traces(3,:)',sum(matAttOnOff)/numCell'
    end
end

trHist=sum(matAttOnOff)/numCell;
plot(trHist);
% grid on

%{
x = [0 reshape(repmat(1:size(s,2)-1,2,1),1,2*size(s,2)-2) size(s,2)];
x = [x fliplr(x)];
y = reshape(repmat(sum(s)/size(s,1),2,1),1,2*size(s,2));
y = [y zeros(1,2*size(s,2))];
% h = patch(x,y,[0.5 0.5 0.5]);  %grey fill
h = patch(x,y,[0 0 0]);  %black fill
set(h,'edgecolor',[0 0 0]);  %black outline
%}
xlim([0 size(region.traces,2)])
set(gca,'xtick',[]);
% title('Network activity','FontWeight','bold')
ylabel('Fraction of ROIs active')


%raster plot of ROI activity----------------------------------------
ax(2)=subplot(3,1,2);
hold on

if ~isempty(stimulusParams)
    for i=1:numel(stimulusParams)
x1=(stimulusParams{i}.frame_indices(1)/stimulusParams{i}.frame_times(1))*stimulusParams{i}.stimulus_times(1);
x2=(stimulusParams{i}.frame_indices(end)/stimulusParams{i}.frame_times(end))*stimulusParams{i}.stimulus_times(end);
x = [x1; x1; x2; x2];
y = [0; length(region.contours)+1; length(region.contours)+1; 0];
h1 = patch(x,y,[0.7 0.7 0.7]);  %grey fill
set(h1,'edgecolor',[0.7 0.7 0.7]);  %gray outline
    end
end

for c = 1:length(region.onsets)
    for d = 1:length(region.onsets{c})
        if strcmp(dur,'true')
            plot(region.onsets{c}(d):region.offsets{c}(d),repmat(c,1,region.offsets{c}(d)-region.onsets{c}(d)+1),'-','color',[0.3 0.3 0.3]);
        else
            plot(region.onsets{c}(d),c,'color',[0.3 0.3 0.3]);
        end
    end
end
ylim([0 length(region.onsets)+1]);
xlim([0 size(region.traces,2)]);
% set(gca,'xtick',[]);
ylabel('ROI no.')
% set(gca,'ytick',(1:fix(length(region.contours)/100))*100)
set(gca,'ydir','reverse')
box on


%bottom raster plot of ROI fluorescence over time---------------------

ax(3)=subplot(3,1,3);
if strcmp(rmBaseline,'true')
m1 = mean(region.traces,1);
m2 = repmat(m1,size(region.traces,1),1);
imagesc(dfoverf(region.traces - m2));
else
imagesc(nt);
end
%set(gca,'xtick',[]);
ylabel('ROI no.')
% set(gca,'clim',[-0.8 0.8])
% set(gca,'ytick',[1 (1:fix(length(region.contours)/100))*100])
set(gca,'ydir','reverse')
set(gca,'xtick',[1 size(nt,2)]);
colorbar('SouthOutside')
secs = fix(size(nt,2)*region.timeres);
str = num2str(secs);
%box on
set(gca,'xticklabel',{'0', str});
xlabel('Time (sec)')
linkaxes(ax,'x');
zoom xon
% figname = 'rastHist'; fname2 = [fnm(1:end-4) figname '.' 'eps']; saveas(gcf,fname2,'epsc');

function rast2matdur = rast2matdur(spk,endpt,sz)
%rast2mat = rast2mat(spk,endpt,dur)
%   converts a rasterplot into a matrix

s = zeros(length(spk),sz);
for c = 1:length(spk)
   for d = 1:length(spk{c})
      s(c,spk{c}(d):endpt{c}(d)) = 1;
   end
end

rast2matdur = s;

function rast2mat = rast2mat(spk,sz)
%rast2mat = rast2mat(spk)
%   converts a rasterplot into a matrix

s = zeros(size(spk,2),sz);
for c = 1:size(spk,2)
   s(c,spk{c}) = 1;
end

rast2mat = s;

function dfoverf = dfoverf(tr)

bases = repmat(mean(tr,2),1,size(tr,2));
dfoverf = (tr-bases)./(bases+eps);
% if mean(f) == 0
%    dfoverf = f / (max(f)-min(f)+eps);
% else
%    a = (f - median(f))/median(f);
%    dfoverf = a ;
% end
