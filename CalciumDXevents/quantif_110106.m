%quantif_110106

%wave durations-- first onset to last onset----
newdur = region.timeres*(region.waveoffsets-region.waveonsets);

%wave ISI------------------------------------------------------------------
tmpres=region.timeres;
d=region.waveonsets;
e=region.waveoffsets;

% if any(d==0)
%     id=find(d==0);
%     d(id)=1;
% end

%     if ~isempty(region.onsets{c})
%         if ismember(c, region.userdata.gdp)
%                 tmpcell = 'gdp';
%         elseif ismember(c,region.userdata.schmutzr)
%                 tmpcell = 'sch';
%         else
%                 tmpcell = 'other';
%         end
%     tmpcellnumb = c;
    
    d1 = [d size(region.traces,2)];
    e1 = [0 e-d];
    ints=diff([0 d1]) - e1;
    if d(1) ~= 1
       ints=ints(2:end); 
    end
    if d(end) ~= size(region.traces,2)
       ints=ints(1:end-1); 
    end
    ints=ints*tmpres;
    
%     celltype = {};
%     cellnumb = {};
%     for z = 1:numel(ints)
%         celltype{z} = tmpcell;
%         cellnumb{z} = num2str(tmpcellnumb);
%     end

%     ints=ints';
%     celltype=celltype';
%     cellnumb=cellnumb';

    %the following few lines concatenates the rowinfo cell array with our numerical results
%     ints1=num2cell(ints);
%     tmp = {};
%     for m = 1:numel(ints);
%       tmp = [tmp; rowinfo;];
%     end
% 
%     ints2 = [tmp cellnumb celltype ints1];
%     intvls=[intvls; ints2;];
%     end

%--------Fraction of ROIs per wave-----------------------------------------
spk = zeros(size(region.traces));
dec = zeros(size(region.traces));
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end
actvFraction=[];
for i=1:numel(region.waveonsets)
    tmp=numel(find(spk(:,region.waveonsets(i):region.waveoffsets(i)) > 0)) / length(region.contours);
    actvFraction=[actvFraction tmp];
end

%-------wave size----------------------------------------------------------
%(1)image dilation (same as 'Maximum' filter in ImageJ of radius 3 on contour movie image array to fill in gaps between ROIs.   Test fr61 and 887.
%(2)make binary image mask for ea fram (pixels > 0)
%(3)use image mask as input to regionprops to calculate area. Will also use mask for weighted centroid calculation (for wave speed, direction calculations)

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m
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

region.waveareapixels={};
for fr=1:size(A1,3)

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
SE = strel('square', 22); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
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
region.waveareapixels{fr}=[STATS.Area];
end

for i=1:numel(region.waveareapixels)
   wavesizemicrons2(i) = (sum(region.waveareapixels{i}) * region.spaceres^2);
end

wavesizemicrons2(1)/((size(A2,1)*region.spaceres)^2) %as percentage of entire area visible in image

% (STATS(1,1).Area * (region.spaceres^2))  %area in um^2.   %could automatically select the ROI with largest area
% (STATS(1,1).Area * (region.spaceres^2))/((size(A2,1)*region.spaceres)^2) %percentage of entire area visible in image



%{
%repeat
A2=imdilate(A2,SE);
figure;imshow(A2)

%lines
SE1=strel('line', 2, 0);  %2 lines are same as rect
SE2=strel('line', 2, 90);
A2=imdilate(A(:,:,fr),[SE1 SE2]);
figure;imshow(A2)


A2=imerode(A2,SE);
figure; imshow(A2)

c=1;
f={};
[szX,szY] = size(region.image);
szZ=size(region.traces,2);
rXY=szY/szX;
szX=szY;
ps = round(region.contours{c});
ps=[ps(:,1) rXY*ps(:,2)];
[subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
inp = inpolygon(subx,suby,region.contours{c}(:,1),rXY*region.contours{c}(:,2));
fx = subx(inp==1);
fy = suby(inp==1);
f{c} = sub2ind([szX, szY],fy,fx);
strel_sz=[numel(unique(fx)) numel(unique(fx))];
%}


%--wave speed---------------------------------------
%(1)Use similar dilation/erosion methods as for wave size calculation and make array of images masks for each wave
%(2)use image masks as input to regionprops to calculate weigthed centroids.

%Otsu's threshold for making wave center mask
% tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
% level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m . Use the with 'includebaseline = true' method. 
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

%--this gives us the square dimensions in pixels of the ROIs (strel_sz)
c=1;
f={};
[szX,szY] = size(region.image);
szZ=size(region.traces,2);
rXY=szY/szX;
szX=szY;
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



%getWaveCentroids2.m
%***Now setup contour image array A againg from myMakeContourMovieWaves.m with 'includebaseline = false' and pass this A array to the wavecentroid script below.

wfr=9;
frames=region.waveonsets(wfr):region.waveoffsets(wfr);
% wavecentroid={}; %use structure if you want to keep more than one centroid for each waveframe
wavecentroid=[]; %if you want to keep just the centroid of the max intensity wave area for each waveframe
% BW_all=zeros(size(A(:,:,1)));
BWmasks={};
for d=1:length(frames)
    %     d=1;
    fr=frames(d);
    %     wavecentroid{d}=[];
    %     fr=883; %testing with diff frames
    % figure;imshow(A1(:,:,fr))
    %use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19px contours.
    % SE = strel('rectangle', [2 2]); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
    SE = strel('rectangle', [2 2]);
    % SE = strel('disk',5);
    % A2=imdilate(A1(:,:,wfr),SE);
    A2=imdilate(A(:,:,fr),SE);
    % figure;imshow(A2)
    % figure; imshow(A(:,:,51))
    
    %try imclose
    SE = strel('square', 22); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
    % SE = strel('disk',5);
    %     A2=imclose(A1(:,:,fr),SE);  %this performs the sequential dilation and erosion operations
    % A2=imdilate(A2,SE);
    A2=imclose(A2,SE);
    % figure;imshow(A2)
    
    % level = graythresh(A2); %otsu's threshold method, including all zero pixels-- will separate into two populations
    % tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
    tmp=A1(:,:,wfr); tmp=tmp(:);  %Base otsu threshold on the total range of  maximum non-zero pixel values within the merged waveframe image
    level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--
    BW = im2bw(A2,level);
    % figure; imshow(BW)
    [L, num] = bwlabel(BW, 8);
    % disp(num)
    RGB = label2rgb(L);
    if num == 0
        %         wavecentroid{d}=[wavecentroid{d}; NaN];
        wavecentroid=[wavecentroid; NaN NaN];
        disp([num2str(fr) 'out first'])
        continue %continues onto next iteration of for loop
    end
    
    STATS = regionprops(L, 'Area');
    areas=[STATS.Area];
    idx=find(areas == max(areas));
    % disp(fr);
    % disp([STATS.Area])
    
    if length(idx) > 1  %in case there are multiple detected areas of same size, like 3x size(ROI), then probably means the detected transients in the wave are too noisy to reliably detect the waveframe center anyways.
        wavecentroid=[wavecentroid; NaN NaN];
        disp([num2str(fr) 'out next'])
        continue
    end
    
    if areas(idx(1)) <= 3*prod(strel_sz+1)  %strel_sz is the the dimensions -1px of the rois used in region.contours for analysis. This gets rid of ones that are not big enough to be part of a wave for detecting the true wave center of mass.
        %         wavecentroid{d}=[wavecentroid{d}; NaN];
        wavecentroid=[wavecentroid; NaN NaN];
        disp([num2str(fr) 'out second'])
        continue %continues onto next iteration of for loop
    end
    
%         figure; imshow(RGB)
    %     hold
    STATS = regionprops(L,A2,'WeightedCentroid');
    % figure; imshow(A2)
    %{
    hold on
    for i=1:(length([STATS.WeightedCentroid])/2)
        plot(STATS(i).WeightedCentroid(1),STATS(i).WeightedCentroid(2),'+r')
        %         wavecentroid{d}=[wavecentroid{d}; STATS(i).WeightedCentroid];
        wavecentroid(d,:)=[wavecentroid; STATS(i).WeightedCentroid];
    end
    hold off
    %}
    centroidlist=[];
    for i=idx
        centroidlist=[centroidlist; STATS(i).WeightedCentroid];
    end
    centroidlist=mean(centroidlist,1);
    
    %----------------------------------------------------------------------
    %algorithm for limiting the moving centroid frames to those in the wavefront
    % CC = bwconncomp(BW,8);
    % numPixels = cellfun(@numel,CC.PixelIdxList);
    % [biggest,idx] = max(numPixels);
    % BW2=zeros(size(BW));
    % BW2(CC.PixelIdxList{idx}) = 1;
    % figure; imshow(BW)
    % figure; imshow(BW2)
    
    % figure; imshow(BW)
    BW2=zeros(size(L));  %former
    [I,J]=find(L == idx(1));
    IND=sub2ind(size(L),I,J);
    BW2(IND)=1; %former
    %     BW_all(IND)=1;  %new
    BW2=im2bw(BW2);
    figure; imshow(BW2)
    
%     BW3 = bwperim(BW2); %former
%     %     BW3 = bwperim(BW_all); %new
%     [row,col]=find(BW3);
%     %     [I,J]=ind2sub(size(BW),find(BW3));
%         figure; imshow(BW3)
    hold on;
    
    %----------------------------------------------------------------------
    %{
    inp=[];
    [szX,szY] = size(region.image);
    szZ=size(region.traces,2);
    rXY=szY/szX;
    szX=szY;
    for c=1:length(region.contours)
        % c=20;
        ps = round(region.contours{c});
        ps=[ps(:,1) rXY*ps(:,2)];
        if rXY > 1
            idx=find(ps(:,2) == min(ps(:,2)));
            ps(idx,2)=min(ps(:,2))-rXY;
        end
        
        centr=centroid(ps);
        K=inpolygon(centr(2),centr(1),I,J);  %find those that are bounded by the original ROI image mask
        if K == 1
            inp=[inp c];
        end
    end
    ons=[];
    offs=[];
    for c=1:numel(inp)
        id=find(region.onsets{inp(c)} >= region.waveonsets(wfr) & region.onsets{inp(c)} <= region.waveoffsets(wfr));
        ons=[ons region.onsets{inp(c)}(id)];
        offs=[offs region.offsets{inp(c)}(id)];
    end
    onset=round(mean(ons));
    offset=round(mean(offs));
    avgTrace=mean(region.traces(inp,frames),1);
    [mx,j]=max(avgTrace);
    %}
    
    
    %----------------------------------------------------------------------
    K=0;
    if ~isempty(BWmasks)
        for i = 1:length(BWmasks)
            %             inp=inpolygon(centroidlist(1),centroidlist(2),BWmasks{i}(:,2),BWmasks{i}(:,1));  %find those that are bounded by the original ROI image mask
            if BWmasks{i}(round(centroidlist(2)),round(centroidlist(1)))
                K=1;
            end
%             if inp == 1
%                 %     disp(num2str(K))
%                 K=1;
%             end
        end
    end
    
    disp(num2str(round(centroidlist)))
    
%     disp(num2str(K))
    %here is final decision point
    %     if fr >= onset && fr<= frames(j)
    %     if fr >= onset && fr<= offset
    if K == 0
        wavecentroid = [wavecentroid; centroidlist];
        plot(wavecentroid(d,1),wavecentroid(d,2),'+r')
        disp([num2str(fr) ' in, K = ' num2str(K)])
        BWmasks{length(BWmasks)+1}=BW2;
    else
        plot(centroidlist(1),centroidlist(2),'ob')
        wavecentroid=[wavecentroid; NaN NaN];
        disp([num2str(fr) 'out, K = ' num2str(K)])
        %         continue
    end
    % figure; plot(avgTrace)
    %----------------------------------------------------------------------
    
    % wavecentroid{fr}=[STATS.WeightedCentroid];
end

ind=1:size(wavecentroid,1);
goodind=ind(find(~isnan(wavecentroid(:,1))));

if length(goodind) > 1
diffind=diff(goodind)';
pxdist=sqrt((abs(diff(wavecentroid(goodind,1)))).^2 + (abs(diff(wavecentroid(goodind,2)))).^2);
consecutivespeeds=(pxdist*region.spaceres)./(diffind*region.timeres)
wavespeed=mean(consecutivespeeds)
% pxdist=pxdist(~isnan(pxdist));
else
    disp('not enough detected center of masses')
end



figure; plot(wavecentroid(:,1),wavecentroid(:,2),'o-');
ylim([0 size(A(:,:,1),1)]);
xlim([0 size(A(:,:,2),2)]);
set(gca,'ydir','reverse'); axis square


K = 0
if BWmasks{3}(258,301)
    K = 1
end



%{
STATS = regionprops(L,'Centroid');
figure; imshow(A2)
hold
for i=1:numel(STATS.Centroid)
plot(STATS(i).Centroid(1),STATS(i).Centroid(2),'+r')
end

STATS = regionprops(L,'Orientation');
STATS(1).Orientation

STATS = regionprops(L,'Eccentricity');
STATS(1).Eccentricity


STATS = regionprops(L, 'Area');
region.waveareapixels{fr}=[STATS.Area];
%}


%--wave direction----------------------------------------------------------
%implemented in getWaveDirections.m on Jan 14, 2011.
%polar plot example for wind speeds
%get angle from arctan (atan2) or cart2pol

% region.wavecentroids
% region.wavespeeds_umpersec
% region.wavecentroids{1}
% theta=atan2()

theta=region.waveprops.wavedirection.theta_radians(~isnan(region.waveprops.wavedirection.theta_radians));
rho=region.waveprops.wavedistance_px(~isnan(region.waveprops.wavedirection.theta_radians));
rho = rho .* region.spaceres;
[x,y] = pol2cart(theta,rho);
figure; compass(x,y)
