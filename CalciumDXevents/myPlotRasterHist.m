function myPlotRasterHist(fnm,region,stimuli,rmBaseline,dur,smoothHist,limitToWaves)
%2011 James B. Ackman
%plot spike rasterplot with overlying histogram
%if dur = 'true', then a onset-duration rasterplot and histogram is generated
%if rmBaseline = 'true', then baseline drift correction is performed
%usage: myPlotRasterHist(fnm,region)
%myPlotRasterHist([],region,[],[],[],'false','false')
%myPlotRasterHist([],region,[],[],[],'true','true')
%myPlotRasterHist([],region,region.stimuli,[],[],'false','false')

%James Ackman, 7/24/2009
%added stimulus parameters, James Ackman, 8/2010
%edited raster plot, James Ackman, 1/5/2011
%May 2011.  added 'smoothHist' argin which can be set to 'true'
%Jun 2011. added limitToWaves argin which can be set to 'true'

if nargin < 7, limitToWaves = 'false'; end
if nargin < 6, smoothHist = 'false'; end
if nargin < 5, dur = 'false'; end
if nargin < 4, rmBaseline = 'false'; end
if nargin < 3, stimuli = []; end

%Network activity plots-----------------------------------------------------
%Left histogram plot of ROI activity-------------------------------------
% nt = filtfilt(fir1(10,0.1,'low'),1,region.traces);
% nt = dfoverf(region.traces);
nt = zeros(size(region.traces));
for c = 1:size(region.traces,1)
    nt(c,:) = dfoverf(region.traces(c,:))*100;
end



figure();
scrsize = get(0,'screensize');
% set(gcf,'Position',[scrsize(3)/2-8.5/11*0.86*scrsize(4)/2 0.07*scrsize(4) 8.5/11*0.86*scrsize(4)*0.5 0.86*scrsize(4)]);
set(gcf,'Position',scrsize);
set(gcf,'color',[1 1 1]);
set(gcf,'PaperType','usletter');
set(gcf,'PaperPositionMode','auto');

ax(1)=subplot(3,1,1);
if strcmp(dur,'true')
    s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
else
    s = rast2mat(region.onsets,size(region.traces,2));
% s = zeros(length(region.onsets),size(region.traces,2));
end

% %need to include the following lines in the above if-else construct so that 'dur' option can be properly utilized.
% if strcmp(limitToWaves,'true')
%    s=zeros(size(region.traces));
%     for c=1:size(region.traces,1)
%         locationIndex = region.location(c);
%         if ~isempty(region.onsets{c})
%             tempOnsets = [];
%             tempOffsets = [];
%             for d=1:length(region.onsets{c})
%                 %            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = region.traces(c,region.onsets{c}(d):region.offsets{c}(d));
%                 %             ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = ntFilt(c,region.onsets{c}(d):region.offsets{c}(d));
%                 ind = find(region.wavedata{locationIndex}.waveonsets <= region.onsets{c}(d));
%                 if ~isempty(ind)
%                     if region.onsets{c}(d) >= region.wavedata{locationIndex}.waveonsets(ind(end)) && region.onsets{c}(d) <= region.wavedata{locationIndex}.waveoffsets(ind(end))
%                         s(c,region.onsets{c}(d)) = 1;
%                         tempOnsets = [tempOnsets region.onsets{c}(d)];
%                         tempOffsets = [tempOffsets region.offsets{c}(d)];
%                     end
%                 end
%             end
%             region.onsets{c} = tempOnsets;
%             region.offsets{c} = tempOffsets;
%         end
%     end 
% end

%need to include the following lines in the above if-else construct so that 'dur' option can be properly utilized.
if strcmp(limitToWaves,'true')
   s=zeros(size(region.traces));
    for c=1:size(region.traces,1)
        locationIndex = region.location(c);
        if ~isempty(region.onsets{c})
            tempOnsets = [];
            tempOffsets = [];
%             for d=1:length(region.onsets{c})
            for d = 1:numel(region.wavedata{locationIndex}.waveonsets);
%                 ind = find(region.wavedata{locationIndex}.waveonsets <= region.onsets{c}(d));
                ind = find(region.onsets{c} >= region.wavedata{locationIndex}.waveonsets(d) & region.onsets{c}<= region.wavedata{locationIndex}.waveoffsets(d));
                if ~isempty(ind)
%                     if region.onsets{c}(d) >= region.wavedata{locationIndex}.waveonsets(ind(end)) && region.onsets{c}(d) <= region.wavedata{locationIndex}.waveoffsets(ind(end))
                        s(c,region.onsets{c}(ind(1))) = 1;
                        tempOnsets = [tempOnsets region.onsets{c}(ind(1))];
                        tempOffsets = [tempOffsets region.offsets{c}(ind(1))];
%                     end
                end
            end
            region.onsets{c} = tempOnsets;
            region.offsets{c} = tempOffsets;
        end
    end 
end

%----
% spk=region.onsets;
% sz=size(region.traces,2);
% s = zeros(size(spk,2),sz);
% for c = 1:size(spk,2)
%    s(c,spk{c}) = 1;
% end
% %----
if strcmp(smoothHist,'true')
    %----alternate hist method, not involving patch for smooth Hist curve---------------
    spk = zeros(size(nt));
    dec = zeros(size(nt));
    for c = 1:size(spk,1)
        spk(c,region.onsets{c}) = 1;
        dec(c,region.offsets{c}) = 1;
    end
    
%     onsets = cell(1,size(spk,1));
%     offsets = cell(1,size(dec,1));
%     for c = 1:size(spk,1)
%         onsets{c} = find(spk(c,:)==1);
%         offsets{c} = find(dec(c,:)==1);
%     end
%     
%     region.onsets=onsets;
%     region.offsets=offsets;
    
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
    plot(trHist,'color',[0.3 0.3 0.3]);
    
else
    x = [0 reshape(repmat(1:size(s,2)-1,2,1),1,2*size(s,2)-2) size(s,2)];
    x = [x fliplr(x)];
    y = reshape(repmat(sum(s)/size(s,1),2,1),1,2*size(s,2));
    y = [y zeros(1,2*size(s,2))];
    % h = patch(x,y,[0.5 0.5 0.5]);  %grey fill
    h = patch(x,y,[0 0 0]);  %black fill
    set(h,'edgecolor',[0 0 0]);  %black outline
end

xlim([0 size(region.traces,2)])
set(gca,'xtick',[]);
% title('Network activity','FontWeight','bold')
ylabel('Fraction of ROIs active')


%raster plot of ROI activity----------------------------------------
ax(2)=subplot(3,1,2);
hold on

if ~isempty(stimuli)
    j = 0;
    k = [];
    hAll = [];
    for numStim = 1:numel(stimuli)
%        if numel(stimuli) > 1
%             mycolors = (hsv(numel(stimuli))) .* 1;
%             mycolors = lines(numel(stimuli)) .* 1;
%             mycolors = cool(numel(stimuli));
            mycolors = [0.8 0.8 1.0; 0.8 1.0 0.8; 1.0 0.8 0.8; 0.6 0.6 1.0; 0.6 1.0 0.6; 1.0 0.6 0.6; 0.4 0.4 1.0; 0.4 1.0 0.4; 1.0 0.4 0.4];
 %       else
 %           mycolors = [0.8 0.8 0.8];
 %       end
        
        for i=1:numel(stimuli{numStim}.stimulusParams)
            x1=(stimuli{numStim}.stimulusParams{i}.frame_indices(1)/stimuli{numStim}.stimulusParams{i}.frame_times(1))*stimuli{numStim}.stimulusParams{i}.stimulus_times(1);
            x2=(stimuli{numStim}.stimulusParams{i}.frame_indices(end)/stimuli{numStim}.stimulusParams{i}.frame_times(end))*stimuli{numStim}.stimulusParams{i}.stimulus_times(end);
            x = [x1; x1; x2; x2];
            y = [0; length(region.contours)+1; length(region.contours)+1; 0];
            h1(i) = patch(x,y,mycolors(numStim,:));  % fill
%             set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
            set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
%             set(h1(i),'FaceAlpha',0.1,'EdgeAlpha',0.1)  %looks great but matlab does not export transparency well
            set(h1(i),'DisplayName',stimuli{numStim}.description{1})
        end
        j = j+i;
        k = [k j];
        hAll = [hAll h1];
%         legend(h1{numStim}(1),'Location','Best')
% opengl('OpenGLBitmapZbufferBug',1); opengl('OpenGLWobbleTesselatorBug',1); opengl('OpenGLLineSmoothingBug',1); opengl('OpenGLEraseModeBug',1);
    end
%     legend([h1{1}(1) h1{2}(1)],'Location','Best')
    legend(hAll(k),'Location','Best')
end

if strcmp(dur,'true')
    for c = 1:size(region.traces,1)
        for d = 1:length(region.onsets{c})
            plot(region.onsets{c}(d):region.offsets{c}(d),repmat(c,1,region.offsets{c}(d)-region.onsets{c}(d)+1),'-','color',[0.3 0.3 0.3]);
        end
    end
    
else
    for c = 1:size(region.traces,1)
        for d = 1:length(region.onsets{c})
            %             plot(region.onsets{c}(d),c,'color',[0.3 0.3 0.3]);
            plot(region.onsets{c}(d),c,'ks','MarkerFaceColor','k','MarkerSize',1);
%             plot(region.onsets{c}(d),c,'o','color',[0.3 0.3 0.3],'MarkerSize',5);
        end
    end 
end
locationMarkers = unique(region.location);
for j = 1:numel(locationMarkers)
	idx = find(region.location == locationMarkers(j));
	text(10,idx(5),region.name{locationMarkers(j)},'FontSize',9,'color',[0.5 0.5 0.5],'HorizontalAlignment','right','Rotation',90);
	if j > 1
	line([0 size(region.traces,2)],[idx(1) idx(1)],'LineStyle','--','LineWidth',0.25,'Color',[0.75 0.75 0.75])
	end
end

ylim([0 length(region.contours)+1]);
xlim([0 size(region.traces,2)]);
% set(gca,'xtick',[]);
ylabel('ROI no.')
% set(gca,'ytick',(1:fix(length(region.contours)/100))*100)
set(gca,'ydir','reverse')
box on
xlabel('Time (frames)')


%bottom raster plot of ROI fluorescence over time---------------------

ax(3)=subplot(3,1,3);

if isfield(region,'tracesFilt')
    imagesc(region.tracesFilt);
else
    if strcmp(rmBaseline,'true')
        m1 = mean(region.traces,1);
        m2 = repmat(m1,size(region.traces,1),1);
        imagesc(dfoverf(region.traces - m2));
    else
        imagesc(nt);
    end
end
%set(gca,'xtick',[]);
ylabel('ROI no.')
% set(gca,'clim',[-0.8 0.8])
% set(gca,'ytick',[1 (1:fix(length(region.contours)/100))*100])
set(gca,'ydir','reverse')
set(gca,'xtick',[1 size(nt,2)]);
% colorbar('SouthOutside')
colorbar('East')
secs = fix(size(nt,2)*region.timeres);
str = num2str(secs);
%box on
set(gca,'xticklabel',{'0', str});
xlabel('Time (sec)')
linkaxes(ax,'x');
pan xon
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
