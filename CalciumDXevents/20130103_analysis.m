%---------------Plot centroids of activations for segmented movie-------------------------
figure; imshow(zeros(size(region.image)));  
num=1;  %placeholder from calciumdx gui code    
handlCoord{num} = [];    
hold on  
% for numcoords = 1:length(region.coords)
for numcoords = 1:3;
    if prod(max(region.coords{numcoords})) ~= prod(size(region.image))
        hCoord = plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
        handlCoord{num} = [handlCoord{num} hCoord];
    end  
end
hold on

for i=1:length(STATS)
	centr = STATS(i).Centroid;
	plot(centr(1),centr(2),'oy','MarkerSize',10)
end




%--------------------------To get MAX intensity based movie (from background subtracted, segmented movie) based on connected components --------------------------------------------
%load('/Volumes/Vega/Users/ackman/Data/2photon/120703i/120703_01_AVG_dummy_outlines_hemis.mat')  
%load('/Volumes/Vega/Users/ackman/Data/2photon/120703i/120703_01_connComponents.mat')  
%fnm = '/Volumes/Vega/Users/ackman/Data/2photon/120703i/120703_01_backgroundSubtract.tif'
tic;
load('120703_01_AVG_dummy_outlines_hemis.mat')  
load('120703_01_connComponents.mat')  
fnm = '120703_01_backgroundSubtract.tif'

[data, series1] = myOpenOMEtiff(fnm);
%load('/Volumes/Vega/Users/ackman/Data/2photon/120703i/120703_01_AVG_dummy_outlines_areas.mat');
%waveONidx=240;
%waveOFFidx=280;
%A = double(series1(:,:,waveONidx:waveOFFidx));
A = double(series1);
clear data series1
%Amean = mean(A,3);
%for i = 1:size(A,3)
%%     A(:,:,i) = (A(:,:,i) - region.image)./region.image;
%    A(:,:,i) = (A(:,:,i) - Amean)./Amean;
%end
%tic
%A2 = wholeBrain_segmentation(A,[],region);
%[A3, CC] = wholeBrain_kmeans(A2,A);
%toc
%tic
%A2 = wholeBrain_segmentation(series1,[],region);
%[A3, CC] = wholeBrain_kmeans(A2,series1);
%toc

%Now make a binary movie array based on the segmented functional signal domains
A3 = zeros(size(A),'uint8');
for i = 1:CC.NumObjects
A3(CC.PixelIdxList{i}) = 1;
end
A3 = logical(A3);

%-------Make movie array with with raw intensity values within the functional ROIs--------
%minValue = abs(min(A(:)));
%A4 = (A + minValue);  %.*2^16;  %scale data to fit within uint16 range [0,2^16), because df/F gives negative values   %commented out, keep double for mat2gray
A4 = A;
A4(A3<1) = 0;  %A3 from wholeBrain_kmeans and wholeBrain_segmentation
A4=mat2gray(A4);

A3double = A;
A3double(A3<1) = 0;  %A3 from wholeBrain_kmeans and wholeBrain_segmentation
A3double(A3>0) = 1;  %<--- Uncomment if you want to get normalized frequency information for the stack instead of 
A3double=mat2gray(A3double);
A3count = sum(A3double,3);

%---------------------sum/max projection of possible domains------------------------------
sm=sum(A4,3);
%imshow(sm,'displayrange',[],'colormap',jet(256)); title('sumproj of raw A4 array')  %imshow requires java to run
imagesc(mat2gray(sm)); title('sumproj of raw A4 array'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])


%***contour plot***
contour(mat2gray(flipud(sm))); title('contour plot of sumproj of raw A4 array');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	


%***contour plot***
contour(mat2gray(flipud(sm)),50); title('contour plot of sumproj of raw A4 array, 50levels');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	


mx=max(sm(:));
A5 = sm./mx;
%imshow(A5,'displayrange',[],'colormap',jet(256)); title('sumproj normalized to max(sumproj)') %imshow requires java to run
imagesc(mat2gray(A5)); title('sumproj norm to max(sumproj)'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])

A5 = max(A4,[],3);
%imshow(A5,'displayrange',[],'colormap',jet(256)); title('maxproj of raw A4 array')  %imshow requires java to run
imagesc(mat2gray(A5)); title('maxproj of raw A4 array'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	


%***contour plot***
contour(mat2gray(flipud(A5))); title('contour plot of maxproj of raw A4 array');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	

%***contour plot***
contour(mat2gray(flipud(A5)),50); title('contour plot of maxproj of raw A4 array, 50levels');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	

%***contour plot***
contour(mat2gray(flipud(A5)),10); title('contour plot of maxproj of raw A4 array, 10levels');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	

%***contour plot***
contour(mat2gray(flipud(A5)),20); title('contour plot of maxproj of raw A4 array, 20levels');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	

%***contour plot***
contour(mat2gray(flipud(A5)),30); title('contour plot of maxproj of raw A4 array, 30levels');
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])	




%----------- now get normalized pixel activation frequency instead
mx=max(sm(:));
A5 = sm./A3count;
%imshow(A5,'displayrange',[],'colormap',jet(256)); title('sumproj normalized to max(sumproj)') %imshow requires java to run
imagesc(mat2gray(A5)); title('sumproj norm to A3count of activations at ea loca'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])

imagesc(mat2gray(A3count)); title('A3count, no. of px activations'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])

mx=max(A3count(:));
A5 = A3count./mx;
%imshow(A5,'displayrange',[],'colormap',jet(256)); title('sumproj normalized to max(sumproj)') %imshow requires java to run
imagesc(mat2gray(A5)); title('A3count of px normalized to max A3count'); colorbar('location','eastoutside')
fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
print('-dpng', [fnm2 '.png'])
print('-depsc', [fnm2 '.eps'])


%{
%---------------------Save movie array as avi---------------------------------------------
%A4 = mat2gray(A4);
%A4 = A4.*2^16;  %scale data to fit within uint16 range [0,2^16)
for fr=1:size(A4,3)
I=A4(:,:,fr);
[I2, map] = gray2ind(I, 256); %figure; imshow(I2,map)
M(fr) = im2frame(I2,map);
end

vidObj = VideoWriter([fnm(1:length(fnm)-4) 'wholeBrain_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])  %VideoWriter not available in MATLAB 2009b
open(vidObj)
for i =1:numel(M)
writeVideo(vidObj,M(i));
end
close(vidObj)
%}

fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS')];
save([fnm2 '.mat'],'A4','CC','A5','A3','sm','mx','-7.3')   %v7.3 file need to save >2GB
toc;