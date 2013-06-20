cellnum = round(get(numslider,'Value'));
delete(th(find(th>0)));
set(cnt,'facecolor',[0 0 0]);
set(cnt(cellnum),'facecolor',1-cl(region.location(cellnum),:));
subplot('position',[0.1 0.15 0.74 0.4]);
hold on
th = [];
if get(df_check,'value') == 0
    if get(halo_check,'value') == 1 & ~isempty(region.halotraces)
        th(1) = plot(region.timeres*(0:size(region.halotraces,2)-1),region.halotraces(cellnum,:),'-k');
        xlim([0 region.timeres*(size(region.halotraces,2)-1)])
    end
    if ~isempty(region.traces)
        th(2) = plot(region.timeres*(0:size(region.traces,2)-1),region.traces(cellnum,:));
        xlim([0 region.timeres*(size(region.traces,2)-1)])
    end
    if ~isempty(region.onsets{cellnum})
        th(3) = plot(region.timeres*(region.onsets{cellnum}-1),region.traces(cellnum,region.onsets{cellnum}),'ro');
    end
    if ~isempty(region.offsets{cellnum})
        th(4) = plot(region.timeres*(region.offsets{cellnum}-1),region.traces(cellnum,region.offsets{cellnum}),'go');
    end
    ylabel('Fluorescence (ADU)');
else
    if get(halo_check,'value') == 1 & ~isempty(region.halotraces)
        th(1) = plot(region.timeres*(0:size(region.halotraces,2)-1),100*(region.halotraces(cellnum,:)-mean(region.halotraces(cellnum,:)))/mean(region.halotraces(cellnum,:)),'-k');
        xlim([0 region.timeres*(size(region.halotraces,2)-1)])
    end
    if ~isempty(region.traces)
        th(2) = plot(region.timeres*(0:size(region.traces,2)-1),100*(region.traces(cellnum,:)-mean(region.traces(cellnum,:)))/mean(region.traces(cellnum,:)));
        xlim([0 region.timeres*(size(region.traces,2)-1)])
    end
    if ~isempty(region.onsets{cellnum})
        th(3) = plot(region.timeres*(region.onsets{cellnum}-1),...
            100*(region.traces(cellnum,region.onsets{cellnum})-mean(region.traces(cellnum,:)))/mean(region.traces(cellnum,:)),'ro');
    end
    if ~isempty(region.offsets{cellnum})
        th(4) = plot(region.timeres*(region.offsets{cellnum}-1),...
            100*(region.traces(cellnum,region.offsets{cellnum})-mean(region.traces(cellnum,:)))/mean(region.traces(cellnum,:)),'go');
    end
    ylabel('Fluorescence (% DF/F)')
end
drawnow;
xlabel('Time (sec)')
title(['Cell #' num2str(cellnum) ' of ' num2str(length(region.contours))]);