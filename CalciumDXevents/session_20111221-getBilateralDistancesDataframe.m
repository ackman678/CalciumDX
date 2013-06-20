%load results structure returned from corrWavesBilateral.m and get info for making dataframe in excel----
%copy and paste relevant info into excel workbook to create dBilateralCorr dataframes for analysis in R-------

%print some info from the results file-----------
fIdx = 2;  %change for results fies that were added to iteratively i nta bastch processin script.
disp(results(fIdx).filename)
if ~isempty(results(fIdx).pairs)
    disp('significant temporal pairs at locations: ')   %Are the hemispheres sign t corr?
    disp(results(fIdx).pairs)
else
    disp('no significant temporal pairs')
end
p_val = results(fIdx).countn./results(fIdx).param.numres;
tmpIdx = find(tril(ones(size(p_val)),-1));
disp('temporal p = ')
disp(num2str(p_val(tmpIdx)'))

disp('Combined Spatial corr counts:')
disp(num2str(results(fIdx).res.DistMetricsCountn))
disp('spatial pvals:')
DistMetricsPvals = results(fIdx).res.DistMetricsCountn./results(fIdx).param.numres;
disp(num2str(DistMetricsPvals))
disp('Mahal Spatial corr counts:')
disp(num2str(results(fIdx).res.MahalmetricsCountn))
disp('Eucl Spatial corr counts:')
disp(num2str(results(fIdx).res.EuclmetricsCountn))

disp('Corresponding wavesDB indices:')
disp(vertcat(results(fIdx).obs.wavesDBind(1:length(results(fIdx).obs.waves(1).waveIND)).pair))
disp('Corresponding file wave numbers:')
disp(vertcat(results(fIdx).obs.wavesIND(:).pair))


%get info on metric values to add to different dataframe in excel for dDistMetrics dataframe analysis in R
disp(num2str((results(fIdx).obs.Mahalmetrics.^0.5)'))
disp('---------------------')
disp(num2str(([results(fIdx).res.Mahalmetrics.Mahalmetrics].^0.5)'))
nanmean(results(fIdx).obs.Mahalmetrics.^0.5)
nanmean([results(fIdx).res.Mahalmetrics.Mahalmetrics].^0.5)

disp(num2str((results(fIdx).obs.Euclmetrics.^0.5)'))
disp('---------------------')
disp(num2str(([results(fIdx).res.Euclmetrics.Euclmetrics].^0.5)'))
nanmean(results(fIdx).obs.Euclmetrics.^0.5)
nanmean([results(fIdx).res.Euclmetrics.Euclmetrics].^0.5)



%fix files (more than one result structure save in each for corrWavesBilateral from 2012-01-04
% fIdx = 2;
% results(fIdx).param.fnmbase
% results = results(fIdx);
% save([results.param.fnmbase '.mat'],'results')




%to get results from corrWavesBilateralSimTcorr.m------------------
sum(sum(results.npairs,2))/sum(sum(results.countn,2))


%----2011-12-24--------Make bilateral wave plot of significantly correlated waves---------------
myColors = jet(15);

load('110323_08_defaultROIs_xcorrn_20111216-233439.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:6],1,'sk--',myColors(1:6,:));
Allcentroids(1).normCentroidsArr = normCentroidsArr;

load('110323_05_defaultROIs_xcorrn_20111216-163701.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2 4],0,'sk--',myColors(7:8,:));
Allcentroids(2).normCentroidsArr = normCentroidsArr;

load('110809_02_70post_xcorrn_20111217-230935.mat')
% [normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2 3],0,'ok--',myColors(9:10,:),-15.76);
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2],0,'ok--',myColors(9,:),-15.76,2,[12:29]);
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[3],0,'ok--',myColors(10,:),-15.76,3,[27:53]);
Allcentroids(3).normCentroidsArr = normCentroidsArr;

load('110809_09_162post_xcorrn_20111219-074839.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:5],0,'ok--',myColors(11:15,:));
Allcentroids(4).normCentroidsArr = normCentroidsArr;

title('biggest area, 50%frames, convexHull centroid, fixfr')
plot2svg






%------2012-01-05----Make bilateral wave plot of significantly correlated waves---------------
myColors = jet(15);

load('110808_03_83post_xcorrn_20120104-012438.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:2],1,'ok--',[],-7,0.1);
Allcentroids(5).normCentroidsArr = normCentroidsArr;

load('110808_08_152_post_37deg_xcorrn_20120103-182431.mat')
% [normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:6],0,'ok--',[]);  %wave 6 gives error
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:5],0,'ok--',[],-7,0.1);
Allcentroids(6).normCentroidsArr = normCentroidsArr;

load('110809_07_135post_xcorrn_20120104-100733.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:7],0,'ok--',[],-10,0.1);
Allcentroids(7).normCentroidsArr = normCentroidsArr;

load('110809_08_148post_xcorrn_20120104-103038.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2],0,'ok--',[],-10,0.1);
Allcentroids(8).normCentroidsArr = normCentroidsArr;
% title('improved-2')

load('110809_02_70post_xcorrn_20111217-230935.mat')
% [normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2 3],0,'ok--',myColors(9:10,:),-15.76);
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2],0,'ok--',[],-15.76,0.1,2,[12:29]);
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[3],0,'ok--',[],-15.76,0.1,3,[27:53]);
Allcentroids(3).normCentroidsArr = normCentroidsArr;

load('110809_09_162post_xcorrn_20111219-074839.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:5],0,'ok--',[],0.1);
Allcentroids(4).normCentroidsArr = normCentroidsArr;

load('110323_05_defaultROIs_xcorrn_20111216-163701.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[2 4],0,'sk--',[],0.1);
Allcentroids(2).normCentroidsArr = normCentroidsArr;

load('110323_08_defaultROIs_xcorrn_20111216-233439.mat')
[normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,[1:6],0,'sk--',[],0.5);
Allcentroids(1).normCentroidsArr = normCentroidsArr;


% title('biggest area, min(20fr or 50%frames), 5fr start, convexHull centroid, fixfr')
% title('biggest area, 50%frames, orig param, convexHull centroid, fixfr')
title('origparam, alpha0.1,6waves alpha0.5')
plot2svg






%---Make pretty rasters for showing matching bilateral ROI location and direction relationships
[s, L, D, sthick, Lthick, Dthick, ntSpk,sOnsets] = myMakeSpatialDirRasters(region,fname);

s2 = rasterReorder(s,0);
D2 = rasterReorder(D,0);
L2 = rasterReorder(L,1);


% figure;
% ax(1)=subplot(2,1,1);
% imagesc(s)
% ax(2)=subplot(2,1,2);
% imagesc(s2)
% xlim([600 2120])
% 
% figure;
% ax(1)=subplot(2,1,1);
% imagesc(D)
% ax(2)=subplot(2,1,2);
% imagesc(D2)
% xlim([600 2120])
% 
% figure;
% ax(1)=subplot(2,1,1);
% imagesc(L)
% ax(2)=subplot(2,1,2);
% imagesc(L2)
% linkaxes(ax,'x'); zoom xon
% xlim([600 2120])


figure;
ax(1) = subplot(3,1,1)
myRasterPlot(sOnsets); axis square
ax(2) = subplot(3,1,2)
imagesc(D2); axis square; colorbar('Location', 'East')
ax(3) = subplot(3,1,3)
imagesc(L2); axis square; colorbar('Location', 'East')
linkaxes(ax,'x'); zoom xon;
xlim([600 2120])







% -- make bigmesh contour ROI outline--
[bigmesh] = makeMatchingROImeshgrid(region,locationIndexPair);
figure;
for i = 1:length(bigmesh.contours)
   cnt1 = patch(bigmesh.contours{i}([1:end 1],1),bigmesh.contours{i}([1:end 1],2),[0.5 0.5 0.5]);
    set(cnt1,'EdgeColor',[0.5 0.5 0.5]);
    set(cnt1,'FaceColor',[1 1 1]);
%     set(cnt1,'FaceAlpha',0.3)  %looks great but matlab does not export transparency well
    set(cnt1,'LineWidth',0.5) 
end
axis equal; axis ij

cnt1 = patch(region.coords{2}([1:end 1],1),region.coords{2}([1:end 1],2),[0 0 0]);
set(cnt1,'EdgeColor',[0 0 0]);
set(cnt1,'FaceColor',[1 1 1]);
    set(cnt1,'FaceAlpha',0)  %looks great but matlab does not export transparency well
set(cnt1,'LineWidth',2)
set(cnt1,'LineStyle','--')

cnt1 = patch(region.coords{3}([1:end 1],1),region.coords{3}([1:end 1],2),[0 0 0]);
set(cnt1,'EdgeColor',[0 0 0]);
set(cnt1,'FaceColor',[1 1 1]);
    set(cnt1,'FaceAlpha',0)  %looks great but matlab does not export transparency well
set(cnt1,'LineWidth',2)
set(cnt1,'LineStyle','--')
