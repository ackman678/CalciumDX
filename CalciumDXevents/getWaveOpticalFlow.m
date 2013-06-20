function getWaveOpticalFlow(region)

[data, series1] = myOpenOMEtiff;
%load test_opticalflow.mat
%save('test_opticalflow.mat')
% waveONidx=1103;
% waveOFFidx=1182;
waveONidx=240;
waveOFFidx=280;

% waveONidx=region.wavedata{2}.waveonsets(6);
% waveOFFidx=region.wavedata{3}.waveoffsets(8);
A = double(series1(:,:,waveONidx:waveOFFidx));
Amean = mean(A,3);
for i = 1:size(A,3)
%     A(:,:,i) = (A(:,:,i) - region.image)./region.image;
    A(:,:,i) = (A(:,:,i) - Amean)./Amean;
end
%Acrop = A(:,250:size(region.image,2),:);
%clear A;
%A = Acrop; clear Acrop
img1 = A(:,:,26);
img2 = A(:,:,27);
img1 = gaussSmooth(img1,6,'same');
img2 = gaussSmooth(img2,6,'same');
figure; imshow(img1,[]); figure; imshow(img2,[])


%----Make Binary Mask Movie of the bihemispheric waveframes--------------
figure;
sz = size(A);
clear M

regionMask1 = poly2mask(region.coords{2}(:,1),region.coords{2}(:,2),sz(1),sz(2));
regionMask2 = poly2mask(region.coords{3}(:,1),region.coords{3}(:,2),sz(1),sz(2));
%	figure; imshow(regionMask1); 	figure; imshow(regionMask2);
bothMasks = regionMask1+regionMask2;
%figure; imshow(bothMasks);
for fr = 1:size(A,3)
    img1 = gaussSmooth(A(:,:,fr),6,'same');
%    img2 = A(:,:,fr);
	level = graythresh(img1);
	BW = im2bw(img1,level);
	[L, num] = bwlabel(BW, 8);
	RGB = label2rgb(L);
	fr1ind = intersect(find(bothMasks),find(BW));
	RGB2 = zeros(size(BW));
	RGB2(fr1ind) = img1(fr1ind);
	imshow(RGB2,[])  %TESTING
	M(fr) = getframe;
end

%for fr = 1:size(A,3)
%	img1 = gaussSmooth(A(:,:,fr),6,'same');
%	level = graythresh(img1);
%	BW = im2bw(img1,level);
%    B(:,:,fr) = (BW);
%end
movie(M,30)

vidObj = VideoWriter(['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])
open(vidObj)
for i =1:numel(M)
writeVideo(vidObj,M(i))
end
close(vidObj)

%-----------------------------------------------------------------------------------------
%Test Optic Flow
	
movie(M2,30)
figure; contour(sqrt(sumVx.^2 + sumVy.^2)); set(gca,'ydir','reverse');

vidObj = VideoWriter(['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])
open(vidObj)
for i =1:numel(M2)
writeVideo(vidObj,M2(i))
end
close(vidObj)

%-----------------------------------------------------------------------------------------
%Test XCorr2d
sz = size(A);
figure;
i=0;
dXY1 = [];
dXY2 = [];
for fr = 1:size(A,3)
i=i+1;
%    img1 = A(:,:,fr-1);
%    img2 = A(:,:,fr);
		img1 = A(:,:,fr);
		img2 = A(:,:,sz(3)+1-fr);
        img1 = gaussSmooth(img1,6,'same');
        img2 = gaussSmooth(img2,6,'same');
		level1 = graythresh(img1);
		level2 = graythresh(img2);
		BW1 = im2bw(img1,level);
		BW2 = im2bw(img2,level);
		
		C = xcorrn(img1,img2,'same');
		imagesc(C)
		M1(i) = getframe;
		[m n] = find(C == max(C(:)));
		dXY1 = [dXY1; m n];

		C = xcorrn(BW1,BW2,'same');
		imagesc(C)
		M2(i) = getframe;
		[m n] = find(C == max(C(:)));
		dXY2 = [dXY2; m n];
end
movie(M2,30)
figure; plot(dXY1(:,2),dXY1(:,1),'o'); ylim([1 sz(1)]); xlim([1 sz(2)])
figure; plot(dXY2(:,2),dXY2(:,1),'o'); ylim([1 sz(1)]); xlim([1 sz(2)])

vidObj = VideoWriter(['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])
open(vidObj)
for i =1:numel(M2)
writeVideo(vidObj,M2(i))
end
close(vidObj)

%-----------------------------------------------------------------------------------------
%Test Xcorr2d between hemispheres
sz = size(A);
figure;
i=0;
dXY1 = [];
dXY2 = [];
clear M1 M2

regionMask1 = poly2mask(region.coords{2}(:,1),region.coords{2}(:,2),sz(1),sz(2));
regionMask2 = poly2mask(region.coords{3}(:,1),region.coords{3}(:,2),sz(1),sz(2));
%	figure; imshow(regionMask1); 	figure; imshow(regionMask2);
cropSZ=region.orientation.value(1);
regionMask1crop=regionMask1(1:cropSZ,:);
regionMask2crop=regionMask2(cropSZ+1:sz(1),:);
regionMask2crop=flipud(regionMask2crop);

for fr = 1:size(A,3)
	i=i+1;
%	fr=26;
	img1 = A(:,:,fr);
	img1 = gaussSmooth(img1,6,'same');
	%------Get average image of label in two hemispheres, crop into two matching oriented images based on midline
	fr1=img1(1:cropSZ,:);
	fr2=img1(cropSZ+1:sz(1),:);
	fr2=flipud(fr2);
	figure; imagesc(fr1)
	figure; imagesc(fr2)
%	
	level = graythresh(img1);
	BW1 = im2bw(fr1,level);
	BW2 = im2bw(fr2,level);
	fr1ind = intersect(find(regionMask1crop),find(BW1));
	fr2ind = intersect(find(regionMask2crop),find(BW2));
	
%	BW1 = gaussSmooth(double(BW1),6,'same')
%	BW2 = gaussSmooth(double(BW2),6,'same')
	fr1bw = zeros(size(fr1));
%	fr1bw(find(BW1)) = fr1(find(BW1));
	fr1bw(fr1ind) = fr1(fr1ind);
	fr2bw = zeros(size(fr2));
%	fr2bw(find(BW2)) = fr2(find(BW2));
	fr2bw(fr2ind) = fr2(fr2ind);
%	figure; imshow(fr1bw,[])
%	figure; imshow(fr2bw,[])
	
	[minDim,loca] = min([size(BW1,1) size(BW2,1)]);
	df = abs(diff([size(fr1,1) size(fr2,1)]));
	if loca == 1
	tempImg = fr2bw(df+1:size(fr2bw,1),:);
	clear fr2bw;
	fr2bw = tempImg; clear tempImg;
	else
	tempImg = fr1bw(df+1:size(fr1bw,1),:);
	clear fr1bw;
	fr1bw = tempImg; clear tempImg;	
	end
	
%	C = normxcorrn(fr1,fr2,'same');
%	imagesc(C)
%	M1(i) = getframe;
%	[m n] = find(C == max(C(:)));
%	dXY1 = [dXY1; m n];

	C = xcorrn(fr1bw,fr2bw,'same');
	imagesc(C)
	M2(i) = getframe;
	[m n] = find(C == max(C(:)));
	disp(fr);
	disp(max(C(:)));
	disp(numel(m));
	if numel(m) > 1
	dXY2 = [dXY2; 1 1];
	else
	dXY2 = [dXY2; m n];
	end
end
movie(M2,30)

%figure; mycolors = jet(size(dXY1,1)); 
%hold on;
%for i = 1:size(dXY1,1)
%plot(dXY1(i,2),dXY1(i,1),'o','Color',mycolors(i,:));
%end 
%ylim([1 size(fr1,1)]); xlim([1 size(fr1,2)]); hold off;

figure; 
mycolors = jet(size(dXY2,1));
hold on;
for i = 1:size(dXY2,1)
plot(dXY2(i,2),dXY2(i,1),'o','Color',mycolors(i,:));
hold on
end
ylim([1 size(fr1,1)]); xlim([1 size(fr1,2)]); hold off;

vidObj = VideoWriter(['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])
open(vidObj)
for i =1:numel(M2)
writeVideo(vidObj,M2(i))
end
close(vidObj)

% %---the following is equal to a Gaussian blur with sigma = 6 in ImageJ
% w=fspecial('gaussian',[51 51],6);
% IF=imfilter(img1,w,'replicate');
% figure; imagesc(IF); title('w=fspecial("gaussian",[51 51],6));')
% figure; imagesc(img1); title('raw data')
%-----------------------------------------------------------------------------------------

D = pdist2(X,Y,'mahalanobis')

function getRawWaveframes

function shuffleWaveISI

function 

function getWaveOpticalFlow(img1,img2)
% img1=A3(:,:,1110); img2=A3(:,:,1120);
figure; imshow(img1,[]); figure; imshow(img2,[])
figure; [Vx,Vy,reliab] = optFlowLk( img1, img2, [], 4, 3, 3e-6, gcf );

figure; [Vx,Vy,reliab] = optFlowLk( img1, img2, [], 12, 1, 3e-6, gcf ); %this one is nice
figure; imagesc(sqrt(Vx.^2 + Vy.^2)); set(gca,'ydir','reverse');
figure; imagesc(atan2(Vy,Vx)); set(gca,'ydir','reverse');
thetaVxy = atan2(Vy,Vx);  %distance
rhoVxy = sqrt(Vx.^2 + Vy.^2);  %magnitude
%figure; imagesc(1 - (mat2gray(thetaVxy)./(rhoVxy.^2)))
%figure; imagesc(sqrt(Vx.^2 + Vy.^2) .* (atan2(Vy,Vx)+(2*pi))); set(gca,'ydir','reverse');