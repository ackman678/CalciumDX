function [normCentroidsArr] = makeBilateralCorrPlot2(results,wavesDBfnm,iPairs,createFig,linestyle,myColors,rotateImage,opacity,fixWave,goodframes)
%makeBilateralCorrPlot(results,'wavesDB.mat',[1:6],1,'sk--')
%region is your region data structure
%results is the results structure returned from corrWavesBilateral.m
%wavesDBfnm is the filename of your wavesDB return from getWaveMasks.m
%iPairs is a list of indices for sign pairs to plot (results from corrWavesBilateral.m, in notebook.key)
%fighandle is the figure handlde from a previous run of makeBilateralCorrPlot.m if you are adding data to the plot.
%2011-12-23 James B. Ackman
if nargin < 10 || isempty(goodframes), goodframes= []; end
if nargin < 9 || isempty(fixWave), fixWave = 0; end
if nargin < 8 || isempty(opacity), opacity = 0.3; end
if nargin < 7 || isempty(rotateImage), rotateImage = 0; end   %amount to rotate image if acquired image was not acquired with midline perfectly at 0degs. Negative values give cw rotate about center axis using imrotate
if nargin < 6 || isempty(myColors), myColors = jet(length(iPairs)); end   %RGB color array
if nargin < 5 || isempty(linestyle), linestyle = 'sk--'; end %optional line styles input
if nargin < 4 || isempty(createFig), createFig = 1; end   %if overplotting multiple waves from multiple 'results' structures
if nargin < 2 || isempty(wavesDBfnm), wavesDBfnm = 'wavesDB.mat'; end

disp(results.param.fnmbase)

try
    matObj = matfile(wavesDBfnm);
catch exception
    rethrow(exception)
end

if createFig > 0
    fig1=figure;
end
j = 0;
% for i = 1:numel(results.obs.wavesDBind)
for i = iPairs
    j = j+1;
    disp(['i = ' num2str(i)])
    pairIND = results.obs.wavesDBind(:,i).pair;
    %     if significant %find spatial corr sign waves
    wave1data=matObj.imagedata(1,pairIND(1));    %fetch waves amask from wavesDB
    wave2data=matObj.imagedata(1,pairIND(2));    %fetch waves amask from wavesDB
    wave1regioncoords=matObj.regioncoords(1,pairIND(1));    %fetch waves amask from wavesDB
    wave2regioncoords=matObj.regioncoords(1,pairIND(2));    %fetch waves amask from wavesDB
    wave1orientation = matObj.orientation(1,pairIND(1));  %get the orientation values from DB for cropping
    wave2orientation = matObj.orientation(1,pairIND(2));
    
    sz1 = size(wave1data.wavemask) %we will set the no. of comparision frames to the length of the shorter wave array
    sz2 = size(wave2data.wavemask) %we will set the no. of comparision frames to the length of the shorter wave array
    szZ = min([sz1(3) sz2(3)]);
    
    szZ = ceil(szZ/2);  %use half of frames (beginning frames) to limit to the spatially matched period of the waves...
%     szZ = min([20 szZ]);
    
    if i == fixWave && ~isempty(goodframes)
        B1 = wave1data.wavemask(:,:,goodframes);  %binary wave mask array for wave1
        clear wave1data
        B2 = wave2data.wavemask(:,:,goodframes);  %binary wave mask array for wave2
        clear wave2data
    else
        B1 = wave1data.wavemask(:,:,1:szZ);  %binary wave mask array for wave1
        clear wave1data
        B2 = wave2data.wavemask(:,:,1:szZ);  %binary wave mask array for wave2
        clear wave2data
    end
    [centroidsXY, meanCentroidXY, normCentroids, normCentroidsAll, normMaskPerim, normCentroidConvexHull] = getCentroid(B1,B2,wave1regioncoords,wave2regioncoords,wave1orientation,wave2orientation,rotateImage);
    
    %     figure(fig1)
    %plot bilateral mean centroids connected by dotted line, make dotted lines different than the midline one
    %     plot([normCentroids(1).centroids(1) normCentroids(2).centroids(1)], [normCentroids(1).centroids(2) normCentroids(2).centroids(2)],linestyle)  %plot circles if dye==OGB1AM, or squares if dye== CaGreenDextran
    
    cn1 = centroid(normMaskPerim(1).centroids);
    cn2 = centroid(normMaskPerim(2).centroids);
    
    plot([cn1(1) cn2(1)], [cn1(2) cn2(2)],linestyle,'LineWidth',1,'MarkerSize',6)  %plot circles if dye==OGB1AM, or squares if dye== CaGreenDextran
    
    xlim([0 1])
    ylim([-1 1])
    axis square
    hold on
    
    
    %     %--contour of iterative centroids outline---------------------------------------------
    %centroid path connected lines
    %     plot(normCentroidsAll(1).centroids(:,1),normCentroidsAll(1).centroids(:,2),'.--','MarkerEdgeColor',myColors(j,:),'MarkerFaceColor',myColors(j,:));
    %     plot(normCentroidsAll(2).centroids(:,1),normCentroidsAll(2).centroids(:,2),'.--','MarkerEdgeColor',myColors(j,:),'MarkerFaceColor',myColors(j,:));
    
    %centroid outline patch objects
    %     cnt1 = patch(normCentroidsAll(1).centroids([1:end],1),normCentroidsAll(1).centroids([1:end],2),myColors(j,:));
    %     set(cnt1,'edgecolor',[0.5 0.5 0.5]);
    %
    %     cnt2 = patch(normCentroidsAll(2).centroids([1:end],1),normCentroidsAll(2).centroids([1:end],2),myColors(j,:));
    %     set(cnt2,'edgecolor',[0.5 0.5 0.5]);
    
    %mask perim outline patch objects
    cnt1 = patch(normMaskPerim(1).centroids([1:end 1],1),normMaskPerim(1).centroids([1:end 1],2),myColors(j,:));
    set(cnt1,'EdgeColor',[0 0 0]);
    set(cnt1,'FaceAlpha',opacity)  %looks great but matlab does not export transparency well
    set(cnt1,'LineWidth',0.5)
    
    cnt2 = patch(normMaskPerim(2).centroids([1:end 1],1),normMaskPerim(2).centroids([1:end 1],2),myColors(j,:));
    set(cnt2,'EdgeColor',[0 0 0]);
    set(cnt2,'FaceAlpha',opacity)  %looks great but matlab does not export transparency well
    set(cnt2,'LineWidth',0.5)
    
    drawnow
    
    
    normCentroidsArr(i).normCentroids = normCentroids;
end

if createFig > 0
    line([0 1],[0 0],'color', [0.5 0.5 0.5]) %dotted line representing midline
end
end %end function

function [centroidsXY, meanCentroidXY, normCentroids, normCentroidsAll, normMaskPerim, normCentroidConvexHull] = getCentroid(B1,B2,wave1regioncoords,wave2regioncoords,wave1orientation,wave2orientation,rotateImage)

%--------for converting to normalized AP [0 1] and L-M-L [1 -1] coords----------------
maxCoords1 = max(wave1regioncoords.regioncoords);
minCoords1 = min(wave1regioncoords.regioncoords);

maxCoords2 = max(wave2regioncoords.regioncoords);
minCoords2 = min(wave2regioncoords.regioncoords);

mnX = min([minCoords1(1) minCoords2(1)]);
mxX = max([maxCoords1(1) maxCoords2(1)]);

minCoords1(1) = mnX;
minCoords2(1) = mnX;
maxCoords1(1) = mxX;
maxCoords2(1) = mxX;

%     mnY = min([minCoords1(2) minCoords2(2)]);
%     mxY = max([maxCoords1(2) maxCoords2(2)]);
%
%     minCoords1(2) = mnY;
%     minCoords2(2) = mnY;
%     maxCoords1(2) = mxY;
%     maxCoords2(2) = mxY;


%     if strcmp(filename,'110809_02_70post_xcorrn_20111217-230935.mat')
%
%         rot 30degs
%     end

%---------get the centroids-------------------------
clear centroidsXY
clear normCentroidsAll
clear normMaskPerim

%**-----First Hemisphere-----**
centroidsXY(1).centroids = [];
normCentroidsAll(1).centroids = [];
normMaskPerim(1).centroids = [];
sumMask1 = zeros(size(B1,1),size(B1,2));
for fr = 1:size(B1,3)
    %get centroid of larrgest conn component in wavemask
    %regionprops
    template=mat2gray(B1(:,:,fr));
    if rotateImage ~= 0  %just for '110809_02_70post_xcorrn_20111217-230935.mat', which was rotated ccw ~15.76degs during acquisition
        I = imrotate(template,rotateImage,'bilinear','crop');
        template = I;
    end
    
    [L, num] = bwlabel(im2bw(template), 8);  %label connected regions
    STATS = regionprops(L, 'Area','Centroid'); %find areas
    if ~isempty(STATS)  %skip blank frames
        areas=[STATS.Area]; %areas as single vector
        idx=find(areas == max(areas));
        idx = idx(1);
        %add centorid xy coords to counter for wave
        centroidsXY(1).centroids = [centroidsXY(1).centroids; STATS(idx).Centroid];
        [xlocanorm, ylocanorm] = getNormalizedXYcoords(STATS(idx).Centroid, maxCoords1, minCoords1,0);
        normCentroidsAll(1).centroids = [normCentroidsAll(1).centroids; [xlocanorm ylocanorm]];
        %             sumMask1 = sumMask1 + template;
        tmp = zeros(size(L));
        tmp(L == idx) = 1;
        sumMask1 = sumMask1 + tmp;
    end
    
    %         %use mean of centroids instead-------------------------------
    % %         STATS = regionprops(L,template,'WeightedCentroid');
    %             STATS = regionprops(L,'Centroid');
    %         if ~isempty(STATS)
    %             wavecentroids = [];
    %             for i = 1:numel(STATS)
    % %                 wavecentroids = [wavecentroids; STATS(i).WeightedCentroid;];
    %                 wavecentroids = [wavecentroids; STATS(i).Centroid;];
    %             end
    % %                         centroidsXY(1).centroids = [centroidsXY(1).centroids; STATS(idx).Centroid];
    %             centroidsXY(1).centroids = [centroidsXY(1).centroids; mean(wavecentroids,1)];
    %         end
end
[L, num] = bwlabel(im2bw(sumMask1), 8);  %label connected regions
STATS = regionprops(L, 'Area','Centroid'); %find areas
areas=[STATS.Area]; %areas as single vector
idx=find(areas == max(areas));   %if not overlapping w last skip
idx = idx(1);
tmp = zeros(size(L));
tmp(L == idx) = 1;

%     figure; imshow(tmp); title('sumMask1')
CH1 = bwconvhull(im2bw(tmp));
%     figure; imshow(CH1); title('CH1')
%     BP1 = bwperim(CH1);/
[BP1,L] = bwboundaries(CH1,'noholes');
%     figure; imshow(BP1); title('Bperim1')
%     [m n] = find(BP1);
boundary = BP1{1};
%     locatmp = [n m];
locatmp = [boundary(:,2) boundary(:,1)];
for fr = 1:size(locatmp,1)
    locapx = locatmp(fr,:);
    [xlocanorm, ylocanorm] = getNormalizedXYcoords(locapx, maxCoords1, minCoords1,0);
    normMaskPerim(1).centroids = [normMaskPerim(1).centroids; [xlocanorm ylocanorm]];
end

[L, num] = bwlabel(CH1, 8);  %label connected regions
STATS = regionprops(L,'Centroid');
[xlocanorm, ylocanorm] = getNormalizedXYcoords([STATS(1).Centroid], maxCoords1, minCoords1,0);
normCentroidConvexHull(1).centroids = [xlocanorm ylocanorm];
%     centroidsXY(1).centroids



%**-----Second Hemisphere-----**
centroidsXY(2).centroids = [];
normCentroidsAll(2).centroids = [];
normMaskPerim(2).centroids = [];
sumMask2 = zeros(size(B2,1),size(B2,2));
for fr = 1:size(B2,3)
    %get centroid of largest conn component in wavemask
    %regionprops
    template=mat2gray(B2(:,:,fr));
    if rotateImage ~= 0 %just for '110809_02_70post_xcorrn_20111217-230935.mat', which was rotated ccw ~15.76degs during acquisition
        I = imrotate(template,rotateImage,'bilinear','crop');
        template = I;
    end
    [L, num] = bwlabel(im2bw(template), 8);  %label connected regions
    STATS = regionprops(L, 'Area','Centroid'); %find areas
    if ~isempty(STATS)  %skip blank frames
        areas=[STATS.Area]; %areas as single vector
        idx=find(areas == max(areas));
        idx = idx(1);
        %add centorid xy coords to counter for wave
        centroidsXY(2).centroids = [centroidsXY(2).centroids; STATS(idx).Centroid];
        [xlocanorm, ylocanorm] = getNormalizedXYcoords(STATS(idx).Centroid, maxCoords2, minCoords2,1);
        normCentroidsAll(2).centroids = [normCentroidsAll(2).centroids; [xlocanorm ylocanorm]];
        %             sumMask2 = sumMask2 + template;
        tmp = zeros(size(L));
        tmp(L == idx) = 1;
        sumMask2 = sumMask2 + tmp;
    end
    
    %         %use mean of centroids instead-------------------------------
    % %         STATS = regionprops(L,template,'WeightedCentroid');
    %         STATS = regionprops(L,'Centroid');
    %         if ~isempty(STATS)
    %             wavecentroids = [];
    %             for i = 1:numel(STATS)
    % %                 wavecentroids = [wavecentroids; STATS(i).WeightedCentroid;];
    %                 wavecentroids = [wavecentroids; STATS(i).Centroid;];
    %             end
    % %                         centroidsXY(1).centroids = [centroidsXY(1).centroids; STATS(idx).Centroid];
    %             centroidsXY(2).centroids = [centroidsXY(2).centroids; mean(wavecentroids,1)];
    %         end
end
[L, num] = bwlabel(im2bw(sumMask2), 8);  %label connected regions
STATS = regionprops(L, 'Area','Centroid'); %find areas
areas=[STATS.Area]; %areas as single vector
idx=find(areas == max(areas));
idx = idx(1);
tmp = zeros(size(L));
tmp(L == idx) = 1;

%     figure; imshow(tmp); title('sumMask2')
CH2 = bwconvhull(im2bw(tmp));
%     figure; imshow(CH2); title('CH2')
%     BP2 = bwperim(CH2);
[BP2,L] = bwboundaries(CH2,'noholes');
%     figure; imshow(BP2); title('Bperim2')
%     [m n] = find(BP2);
boundary = BP2{1};
%     locatmp = [n m];
locatmp = [boundary(:,2) boundary(:,1)];
for fr = 1:size(locatmp,1)
    locapx = locatmp(fr,:);
    [xlocanorm, ylocanorm] = getNormalizedXYcoords(locapx, maxCoords1, minCoords1,1);
    normMaskPerim(2).centroids = [normMaskPerim(2).centroids; [xlocanorm ylocanorm]];
end

[L, num] = bwlabel(CH2, 8);  %label connected regions
STATS = regionprops(L,'Centroid');
[xlocanorm, ylocanorm] = getNormalizedXYcoords([STATS(1).Centroid], maxCoords2, minCoords2,1);
normCentroidConvexHull(2).centroids = [xlocanorm ylocanorm];
%     centroidsXY(2).centroids







%find mean centroid xy) loca adn add to counter
meanCentroidXY(1).centroid = mean(centroidsXY(1).centroids,1);
meanCentroidXY(2).centroid = mean(centroidsXY(2).centroids,1);

%--------convert mean centroid to normalized AP [0 1] and L-M-L [1 -1] coords----------------
locapx = meanCentroidXY(1).centroid;
[xlocanorm, ylocanorm] = getNormalizedXYcoords(locapx, maxCoords1, minCoords1,0);
normCentroids(1).centroids = [xlocanorm ylocanorm];

locapx = meanCentroidXY(2).centroid;
[xlocanorm, ylocanorm] = getNormalizedXYcoords(locapx, maxCoords2, minCoords2,1);
normCentroids(2).centroids = [xlocanorm ylocanorm];

normCentroids(1).centroids
normCentroids(2).centroids


end %end function

function [xlocanorm, ylocanorm] = getNormalizedXYcoords(locapx, maxCoords, minCoords,flipcoords)
if nargin < 4 || isempty(flipcoords), flipcoords = 0; end
xlocapx = locapx(1);
ylocapx = locapx(2);
if flipcoords < 1
    %     if wave1orientation.orientation.value(1) < mn(2)
    %     xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
    %     ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
    %     ylocanorm = -1*abs(1-ylocanorm);
    %     else
    
    xlocanorm = (xlocapx - minCoords(1))/(maxCoords(1) - minCoords(1));
    ylocanorm = (ylocapx - minCoords(2))/(maxCoords(2) - minCoords(2));
    ylocanorm = abs(1-ylocanorm);
    %     end
    
else
    %     mx = max(wave2regioncoords.regioncoords);
    %     mn = min(wave2regioncoords.regioncoords);
    %     if wave2orientation.orientation.value(1) < mn(2)
    %     else
    %     xlocanorm = (xlocapx - mn(1))/(mx(1) - mn(1));
    %     ylocanorm = (ylocapx - mn(2))/(mx(2) - mn(2));
    %     end
    
    xlocanorm = (xlocapx - minCoords(1))/(maxCoords(1) - minCoords(1));
    ylocanorm = (ylocapx - minCoords(2))/(maxCoords(2) - minCoords(2));
    ylocanorm = -1*abs(1-ylocanorm);
end
end