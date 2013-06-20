%getWaveDistance.m
%James Ackman, 1/17/2011
%{
%--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
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

%}

region.waveprops.eccentricity=[];
region.waveprops.waveorientation_radians=[];
region.waveprops.wavedistance_px=[];

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
    
    CC = bwconncomp(BW,8);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    BW2=zeros(size(BW));
    BW2(CC.PixelIdxList{idx}) = 1;
    %     figure; imshow(BW)
%     figure; imshow(BW2)
    CC = bwconncomp(BW2,8);
    
    % [L, num] = bwlabel(BW, 8);
    % % disp(num)
    % RGB = label2rgb(L);
    % figure; imshow(RGB)
    %
    % STATS = regionprops(L, 'Area');
    % region.waveareapixels{fr}=[STATS.Area];
    
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
    
    region.waveprops.eccentricity=[region.waveprops.eccentricity waveeccent];
    region.waveprops.wavedistance_px=[region.waveprops.wavedistance_px wavedist];
    region.waveprops.waveorientation=[region.waveprops.waveorientation_radians waveori];
    
%     [x,y] = pol2cart(waveori,wavedist);
%     figure; compass(x,y)
%     title(['orientation wave no. ' num2str(fr)])
    
    %{
    %compare with distance from the weighted centroids used for speed and direction calculations. Ellipse diameter calculation above more faithfull represents the wave distance.
    wfr=fr;
    wavecentroid=region.wavecentroids{wfr};  %save data in region structure array
    
    %printout some results. If not enough detected waveframe centroids, no wavespeeds will be calculated.
    ind=1:size(wavecentroid,1);
    goodind=ind(find(~isnan(wavecentroid(:,1))));
    
    if length(goodind) > 1
        diffind=diff(goodind)';
        pxdist=sqrt((abs(diff(wavecentroid(goodind,1)))).^2 + (abs(diff(wavecentroid(goodind,2)))).^2);
%         consecutivespeeds=(pxdist*region.spaceres)./(diffind*region.timeres)
        wavespeed=mean(consecutivespeeds);
        %             disp(['mean wavespeed: ' num2str(wavespeed) ' um/sec'])
        % pxdist=pxdist(~isnan(pxdist));
        %             region.wavespeeds_umpersec{wfr}=consecutivespeeds;
%         disp(num2str(pxdist))
        pxdist_microns=pxdist.*region.spaceres;
%         disp(num2str(pxdist_microns))
        pxdist_microns=sum(pxdist.*region.spaceres);
        disp(['sum dist. ' num2str(pxdist_microns)])
        disp(' ')
    else
        disp('not enough detected center of masses')
        disp(' ')
        %             region.wavespeeds_umpersec{wfr}=NaN;
    end
    %}
    
end
