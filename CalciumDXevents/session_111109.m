%session 2011-11-09
%testing some new image processing techniques for similarity/corr measure between wave patterns

%test example image data loaded from corrspatialtest_111108.mat
load corrspatialtest_111108.mat
BW = im2bw(A3(:,:,1114));  %a bilateral wave frame that has been smoothed and dilated from myMakeContourMoviesWaves.m
figure; imshow(BW); title('A3(:,:,1114)')
D = bwdist(~BW);
figure, imshow(D,[]), title('Distance transform of bw')  %Euclidean

%{
%Test example from documentation, two overlapping circles-----------------
center1 = -10; 
center2 = -center1; 
dist = sqrt(2*(2*center1)^2); 
radius = dist/2 * 1.4; 
lims = [floor(center1-1.2*radius) ceil(center2+1.2*radius)]; 
[x,y] = meshgrid(lims(1):lims(2)); 
bw1 = sqrt((x-center1).^2 + (y-center1).^2) <= radius; 
bw2 = sqrt((x-center2).^2 + (y-center2).^2) <= radius; 
bw = bw1 | bw2; 
figure, imshow(bw), title('bw')

D = bwdist(~bw); 
figure, imshow(D,[]), title('Distance transform of ~bw')
%}

%test normxcorr, xcorr on autocorr, subsequent, and dissimilar wave frames
figure; 
subplot(3,1,1)
B=normxcorr2(A3(:,:,1114),A3(:,:,1114));  %normalized xcorr2. corr values between [-1 1].
imagesc(B); title('normxcorr2 same'); axis square; colorbar
subplot(3,1,2)
B=xcorrn(double(A3(:,:,1114)),double(A3(:,:,1114)));  %xcorrn from piotr ImageVideoProcessingToolbox.
imagesc(B); title('xcorrn same'); axis square; colorbar
subplot(3,1,3)
hist(reshape(B,1,numel(B))); title('xcorrn same')


figure; 
subplot(3,1,1)
B=normxcorr2(A3(:,:,1114),A3(:,:,1113));  %normalized xcorr2. corr values between [-1 1]. 
imagesc(B); title('normxcorr2 similar'); axis square; colorbar
subplot(3,1,2)
B=xcorrn(double(A3(:,:,1114)),double(A3(:,:,1113)));  %xcorrn from piotr ImageVideoProcessingToolbox. 
imagesc(B); title('xcorrn similar'); axis square; colorbar
subplot(3,1,3)
hist(reshape(B,1,numel(B))); title('xcorrn similar')


figure; 
subplot(3,1,1)
B=normxcorr2(A3(:,:,1114),A3(:,:,1153));  %normalized xcorr2. corr values between [-1 1].
imagesc(B); title('normxcorr2 dissimilar'); axis square; colorbar
subplot(3,1,2)
% B=convn(double(A3(:,:,1114)),rot90(double(A3(:,:,1153)),2));  %xcorrn from piotr ImageVideoProcessingToolbox. corr values between [-1 1]. This is an autocorr. Many values above 0.5, and close to or at 1. Highest values directly in image center (256,512). Rest of values offset to upperleft quadrant, different than xcorr2
B=xcorrn(double(A3(:,:,1114)),double(A3(:,:,1153)));  %xcorrn from piotr ImageVideoProcessingToolbox. 
imagesc(B); title('xcorrn dissimilar'); axis square; colorbar
subplot(3,1,3)
hist(reshape(B,1,numel(B)),30); title('xcorrn dissimilar')
%--------------------------------------------------------------------------

%------Get average image of label in two hemispheres, crop into two matching oriented images based on midline
sz=size(region.image);
cropSZ=region.orientation.value(1);
fr1=zeros(cropSZ,sz(2));
fr2=zeros(cropSZ+1,sz(2));

fr1ind=sub2ind(sz,1:cropSZ,1:sz(2));
fr2ind=sub2ind(sz,cropSZ+1:sz(1),1:sz(2));
fr1=region.image(1:size(fr1,1),:);
fr2=region.image(cropSZ+1:sz(1),:);
fr2=flipud(fr2);
fr2=fr2(1:size(fr1),:);
figure; imagesc(fr1)
figure; imagesc(fr2)


% fr1 = I1;
% fr2= I2;

%------Now take the two images and test normxcorr on them vs autocorr-->
%the images are similar enough that the corr matrix between them is pretty much the same as the autocorr matrix
figure; 
subplot(3,2,1)
B=normxcorr2(fr1,fr2);  %normalized xcorr2. corr values between [-1 1]. 
imagesc(B); title('normxcorr2 similar'); axis square; colorbar
subplot(3,2,2)
B=normxcorr2(fr1,fr1);  %normalized xcorr2. corr values between [-1 1]. 
imagesc(B); title('normxcorr2 autocorr'); axis square; colorbar

subplot(3,2,3)
B=xcorrn(double(fr1),double(fr2));  %xcorrn from piotr ImageVideoProcessingToolbox. 
imagesc(B); title('xcorrn similar'); axis square; colorbar
subplot(3,2,4)
B=xcorrn(double(fr1),double(fr1));  %xcorrn from piotr ImageVideoProcessingToolbox. 
imagesc(B); title('xcorrn autocorr'); axis square; colorbar

subplot(3,2,5:6)
hist(reshape(B,1,numel(B))); title('xcorrn similar')
%--------------------------------------------------------------------------

%---Test optic flow--------------------------------------------------------
% create square + translated square (B) + rotated square (C)
tst1=zeros(50,50); tst1(16:35,16:35)=1;
tst2=zeros(50,50); tst2(17:36,17:36)=1;
C=imrotate(tst1,5,'bil','crop');
figure; imshow(tst1,[])
figure; imshow(tst2,[])
figure; imshow(C,[])
optFlowLk( tst1, tst2, [], 2, 2, 3e-6, 1 );
optFlowLk( tst1, C, [], 2, 2, 3e-6, 2 );
% compare on stored real images (of mice)
load optFlowData;
img1=I4A; img2=I4B;
figure; imshow(img1,[]); figure; imshow(img2,[])
[Vx,Vy,reliab] = optFlowLk( img1, img2, [], 4, 1.2, 3e-6, 1 );
[Vx,Vy,reliab] = optFlowCorr( img1, img2, 3, 5, 1.2, .01, 2 );
[Vx,Vy] = optFlowHorn( img1, img2, 2, 3 );

% figure(show); clf; im( I1 );
% hold('on'); quiver( Vx, Vy, 0,'-b' ); hold('off');

img1=A3(:,:,1110); img2=A3(:,:,1120);
figure; imshow(img1,[]); figure; imshow(img2,[])
figure; [Vx,Vy,reliab] = optFlowLk( img1, img2, [], 2, 1.2, 3e-6, gcf );



%----test vol visualization
% 
% figure;
% colordef(gcf,'white')
% load mri;
% D = squeeze(D);
% [x y z D] = subvolume(D, [nan nan nan nan nan 4]);
% p = patch(isosurface(x,y,z,D, 5), 'FaceColor', 'red', 'EdgeColor', 'none');
% p2 = patch(isocaps(x,y,z,D, 5), 'FaceColor', 'interp', 'EdgeColor', 'none');
% isonormals(x,y,z,D,p);
% view(3);
% daspect([1 1 .4])
% colormap(gray(100))
% camva(9);
% box on
% camlight(40, 40);
% camlight(-20,-10);
% lighting gouraud