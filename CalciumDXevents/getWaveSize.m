%getWaveSize
%James Ackman, 1/7/2011
%-------wave size----------------------------------------------------------
%(1)image dilation (same as 'Maximum' filter in ImageJ of radius 3 on contour movie image array to fill in gaps between ROIs.   Test fr61 and 887.
%(2)make binary image mask for ea fram (pixels > 0)
%(3)use image mask as input to regionprops to calculate area. Will also use mask for weighted centroid calculation (for wave speed, direction calculations)

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m

A = myMakeContourMovieWaves([],region,'false');

%-----Now setup a small array containing a single merged frame representing each wave

A1=zeros([size(A(:,:,1)) numel(region.waveonsets)]);
for i=1:numel(region.waveonsets)
    Adummy=A(:,:,region.waveonsets(i));
    for d=region.waveonsets(i):region.waveoffsets(i)
%         figure; imshow(A(:,:,d))
       if d ~= region.waveoffsets(i)
           Adummy=max(Adummy,A(:,:,d+1));
       end
    end
%     figure; imshow(Adummy)
    %Add dummy frame to master wave list length of region.waveonsets
    A1(:,:,i)=Adummy;
end

region.waveprops.waveareapixels={};
for fr=1:length(region.waveonsets)

% fr=5; %testing with diff frames
% figure;imshow(A1(:,:,fr))
%use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19px contours.
% SE = strel('rectangle', [2 2]); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
SE = strel('rectangle', [2 2]);
% SE = strel('disk',5);
A2=imdilate(A1(:,:,fr),SE);
% figure;imshow(A2)

% SE = strel('diamond', 22);
% A2=imclose(A2,SE);
% % A2=imclose(A1(:,:,fr),SE);
% figure; imshow(A2)

% %imerode
% SE = strel('rectangle', [23 23]);
% A2=imerode(A2,SE);
% figure; imshow(A2)

% SE=strel('line', 30, 45);
% A2=imclose(A2,SE);
% figure; imshow(A2)

%try imclose
SE = strel('square', strel_sz(1)+3); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
% SE = strel('disk',5);
% A2=imclose(A1(:,:,fr),SE);  %this performs the sequential dilation and erosion operations
A2=imclose(A2,SE);
% figure;imshow(A2)

level = graythresh(A2); %otsu's threshold method, including all zero pixels-- will separate into two populations
% tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
% level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--
BW = im2bw(A2,level);
% figure; imshow(BW)
[L, num] = bwlabel(BW, 8);
% disp(num)
RGB = label2rgb(L);
figure; imshow(RGB)

STATS = regionprops(L, 'Area');
region.waveprops.waveareapixels{fr}=[STATS.Area];
end

for i=1:numel(region.waveprops.waveareapixels)
   wavesizemicrons2(i) = (sum(region.waveprops.waveareapixels{i}) * region.spaceres^2);
end
disp(num2str(wavesizemicrons2))

% wavesizemicrons2(1)/((size(A2,1)*region.spaceres)^2) %as percentage of entire area visible in image

% (STATS(1,1).Area * (region.spaceres^2))  %area in um^2.   %could automatically select the ROI with largest area
% (STATS(1,1).Area * (region.spaceres^2))/((size(A2,1)*region.spaceres)^2) %percentage of entire area visible in image