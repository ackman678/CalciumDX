function myPlotMultiTrace(region,stimulusParams,spk,dec)
%myPlotMultiTrace
if nargin < 2, stimulusParams = []; end
if nargin < 3, spk = []; end
if nargin < 4, dec = []; end
% nt = dfoverf(region.traces);
tr=region.traces;
bases = repmat(mean(tr,2),1,size(tr,2));
nt = (tr-bases)./(bases+eps);
nt=nt.*100;

figure();
scrsize = get(0,'screensize');
set(gcf,'position',[scrsize(3)/2-8.5/11*0.86*scrsize(4)/2 0.07*scrsize(4) 8.5/11*0.86*scrsize(4)*0.5 0.86*scrsize(4)]); %switch the last two elements in the vector to make landscape orientation
set(gcf,'color',[1 1 1]);
set(gcf,'PaperType','usletter');
set(gcf,'PaperPositionMode','auto');

%To find active cells
% s = rast2mat(region.onsets,size(region.traces,2));
%--------------------------------------------------------------------------
hold on
%calculate spacing in between ea trace to fit all on one plot
dy = abs(max(nt,[],2)) + abs(min(nt,[],2));
dy(1)=0;

%plot stimulus periods
if ~isempty(stimulusParams)
    for i=1:numel(stimulusParams)
        x1=(stimulusParams{i}.frame_indices(1)/stimulusParams{i}.frame_times(1))*stimulusParams{i}.stimulus_times(1);
        x2=(stimulusParams{i}.frame_indices(end)/stimulusParams{i}.frame_times(end))*stimulusParams{i}.stimulus_times(end);
        x = [x1; x1; x2; x2];
        % y = [0; length(region.contours)+1; length(region.contours)+1; 0];
        y = [-20; sum(dy); sum(dy); -20];
        h1 = patch(x,y,[0.7 0.7 0.7]);  %grey fill
        set(h1,'edgecolor',[0.7 0.7 0.7]);  %gray outline
    end
end

%--------------------------------------------------------------------------
%Sort ROIs in ascending order by approximately the time in which the pixels were scanned by the laser.
% function centroid = calciumdxCentroid(coords)
centroid=zeros(numel(region.contours),2);
%{
for i=1:numel(region.contours)
coords=region.contours{i};

if prod(size(coords))==0
   cx = NaN;
   cy = NaN;
else
   m = [coords; coords(1,:)];
   x = m(:,1);
   y = m(:,2);
   
   a = (sum(x(1:end-1).*y(2:end)) - sum(x(2:end).*y(1:end-1)))/2;
   cx = sum((x(1:end-1)+x(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*a);
   cy = sum((y(1:end-1)+y(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*a);
end

centroid(i,:) = [cx cy];
end

centroid(:,2)=centroid(:,2).^2;
[B idx]=sortrows(centroid,2);
%}
idx=1:size(nt,1);
%--------------------------------------------------------------------------
%plot each trace
for j=1:numel(idx)
    i=idx(j);
    %         dy=40;
    nt=nt+dy(i);
    %     plot(nt(i,:),'color',[.75 .75 .75]); %default with gray trace
    plot(nt(i,:),'color',[0.5 0.5 0.5]); %default with gray trace
    
    %plot the event onsets and offsets
    num=i;
    
    if ~isempty(spk)
        f = find(spk(num,:)==1);
        %     cmenu = zeros(1,length(f));
        %     coffmenu = zeros(1,length(f));
        for c = 1:length(f)
            %         cmenu(c) = uicontextmenu;
            %         h6 = uimenu(cmenu(c), 'Label', 'Delete event', 'Callback', ['selev = ' num2str(c) '; hevDeleteEvent;']);
            %         h7 = uimenu(cmenu(c), 'Label', 'Move onset', 'Callback', ['selev = ' num2str(c) '; hevMoveEvent;']);
            %         if c > 1
            %             uimenu(cmenu(c), 'Label', 'Combine with previous', 'Callback', ['selev = ' num2str(c) '; hevCombineEvent;']);
            %         end
            %         h3 = plot(f(c),nt(num,f(c)),'or','uicontextmenu',cmenu(c));
            plot(f(c),nt(num,f(c)),'or');
            
            %         coffmenu(c) = uicontextmenu;
            %         h4 = uimenu(coffmenu(c), 'Label', 'Move offset', 'Callback', ['selev = ' num2str(c) '; hevMoveEventOff;']);
            %         h5 = plot(g(c),nt(num,g(c)),'og','uicontextmenu',coffmenu(c));
        end
    end
    if ~isempty(dec)
        g = find(dec(num,:)==1);
        for c = 1:length(g)
            plot(g(c),nt(num,g(c)),'og');
        end
    end
    
end

ylim([-20 max(nt(i,:))]);
% ylim([-20 sum(dy)]);
xlim([0 size(region.traces,2)]);
% set(gca,'xtick',[]);
% ylabel('Cell no.')
% box on
% linkaxes(ax,'x');
% set(gca,'UserData',[0 size(region.traces,2)])
% set(gca,'buttondownfcn','hevZoom2')

hold off

% for c = 1:length(region.contours)
%     for d = 1:length(region.onsets{c})
%         if strcmp(dur,'true')
%             plot(region.onsets{c}(d):region.offsets{c}(d),repmat(c,1,region.offsets{c}(d)-region.onsets{c}(d)+1),'-','color',[0.3 0.3 0.3]);
%         else
%             plot(region.onsets{c}(d),c,'color',[0.3 0.3 0.3]);
%         end
%     end
% end

% function dfoverf = dfoverf(tr)
% 
% bases = repmat(mean(tr,2),1,size(tr,2));
% dfoverf = (tr-bases)./(bases+eps);