function region = getWaveSizeDistance(region)
%getWaveSizeDistance.m
%James Ackman, 1/17/2011
%updated to work by location, 5/8/2011
%-------wave size----------------------------------------------------------
%(1)image dilation (same as 'Maximum' filter in ImageJ of radius 3 on contour movie image array to fill in gaps between ROIs.   Test fr61 and 887.
%(2)make binary image mask for ea fram (pixels > 0)
%(3)use image mask as input to regionprops to calculate area. Will also use mask for weighted centroid calculation (for wave speed, direction calculations)

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m

locationMarkers = unique(region.location);
regionsAll = splitRegion(region);
for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    results = getWaveSizeDistanceByLocation(tmpregion,locationIndex);
    region.wavedata{locationIndex}.waveprops.waveareapixels=results.waveareapixels;
    region.wavedata{locationIndex}.waveprops.eccentricity=results.eccentricity;
    region.wavedata{locationIndex}.waveprops.waveorientation_radians=results.waveorientation_radians;
    region.wavedata{locationIndex}.waveprops.wavedistance_px=results.wavedistance_px;
end


function results = getWaveSizeDistanceByLocation(region,locationIndex)
% A = myMakeContourMovieWaves([],region,'false');
A = myMakeContourMovieWavesOnsets([],region,locationIndex,'false');

%--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
c=1;
f={};
[szX,szY] = size(region.image);  %assuming szY is the largest dimension
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;  %assuming other rectanqular images (like from CCD camera) don't have evenly divisible number of row lines as column lines (unlike for half scannin with a laser scanning microscope)
end

ps = round(region.contours{c});
ps=[ps(:,1) rXY*ps(:,2)];
if rXY > 1
    idx=find(ps(:,2) == min(ps(:,2)));
    ps(idx,2)=min(ps(:,2))-rXY;
end
[subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
fx = subx(inp==1);
fy = suby(inp==1);
f{c} = sub2ind([szX, szY],fy,fx);
strel_sz=[numel(unique(fx)) numel(unique(fx))];

%-----Now setup a small array containing a single merged frame representing each wave

A1=zeros([size(A(:,:,1)) numel(region.wavedata{locationIndex}.waveonsets)]);
for i=1:numel(region.wavedata{locationIndex}.waveonsets)
    Adummy=A(:,:,region.wavedata{locationIndex}.waveonsets(i));
    for d=region.wavedata{locationIndex}.waveonsets(i):region.wavedata{locationIndex}.waveoffsets(i)
        %         figure; imshow(A(:,:,d))
        if d ~= region.wavedata{locationIndex}.waveoffsets(i)
            Adummy=max(Adummy,A(:,:,d+1));
        end
    end
    %     figure; imshow(Adummy)
    %Add dummy frame to master wave list length of region.wavedata{locationIndex}.waveonsets
    A1(:,:,i)=Adummy;
end

results.waveareapixels={};
results.eccentricity=[];
results.waveorientation_radians=[];
results.wavedistance_px=[];
wavesizemicrons2=[];

for fr=1:length(region.wavedata{locationIndex}.waveonsets)
    
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
%     figure; imshow(RGB)   %FOR TESTING
    
    STATS = regionprops(L, 'Area');
    results.waveareapixels{fr}=[STATS.Area];
    
    CC = bwconncomp(BW,8);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    BW2=zeros(size(BW));
    BW2(CC.PixelIdxList{idx}) = 1;
    %     figure; imshow(BW)
    %     figure; imshow(BW2)
    CC = bwconncomp(BW2,8);
    
    STATS = regionprops(CC, 'all');
    disp(['wave no. ' num2str(fr)])
    %     disp(['area: ' num2str([STATS.Area])]);
    disp(['Eccentricity: ' num2str([STATS.Eccentricity])]);
    disp(['MajorAxisLength: ' num2str([STATS.MajorAxisLength])]);
    diameter_microns= [STATS.MajorAxisLength] .* region.spaceres;
    disp(['MajorAxisLength, microns: ' num2str(diameter_microns)]);
    disp(['Orientation: ' num2str([STATS.Orientation])]);
    
    waveeccent=[STATS.Eccentricity];
    wavedist=[STATS.MajorAxisLength];
    waveori=[STATS.Orientation];
    waveori=waveori.*(pi/180);
    
    results.eccentricity=[results.eccentricity waveeccent];
    results.wavedistance_px=[results.wavedistance_px wavedist];
    results.waveorientation_radians=[results.waveorientation_radians waveori];
end

for i=1:numel(results.waveareapixels)
    wavesizemicrons2(i) = (sum(results.waveareapixels{i}) * region.spaceres^2);
end
disp(num2str(wavesizemicrons2))


% wavesizemicrons2(1)/((size(A2,1)*region.spaceres)^2) %as percentage of entire area visible in image

% (STATS(1,1).Area * (region.spaceres^2))  %area in um^2.   %could automatically select the ROI with largest area
% (STATS(1,1).Area * (region.spaceres^2))/((size(A2,1)*region.spaceres)^2) %percentage of entire area visible in image