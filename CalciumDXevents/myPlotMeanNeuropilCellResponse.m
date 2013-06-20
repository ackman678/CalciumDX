%myPlotMeanNeuropilCellResponse.m
% locationMarkers = unique(region.location);
idx = find(region.location == 1);
trNeuropil = nt(idx,:);
meanNeuropil = mean(trNeuropil,1);
idx = find(region.location == 2);
trCells = nt(idx,:);
meanCells = mean(trCells,1);
figure; 
ax(1) = subplot(4,1,1)
plot(meanNeuropil); ylabel('meanNeuropil, dF/F')
ax(2) = subplot(4,1,2)
imagesc(trNeuropil); ylabel('Roi no.')
ax(3) = subplot(4,1,3)
plot(meanCells); ylabel('meanCells, dF/F')
ax(4) = subplot(4,1,4)
imagesc(trCells); ylabel('Cell no.')
zoom xon
ylim(ax(1),[min(meanCells(:)) max(meanCells(:))]);
ylim(ax(3),[min(meanCells(:)) max(meanCells(:))]);
xlabel('Time (frames)')
