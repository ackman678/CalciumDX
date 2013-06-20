%getWaveCentroids.m
%Used for calculating center of masses in sequential frames during a wave
%The resulting centroid list can be used in wave speed and wave direction calculations----

%revised--2/2/2011  Changed line 92 back to a gray threshold that simply excludes all non-zero pixels from the merged waveframe, instead of threshold based on intensites of all pixels in the image. This is because some images may have bright-non signal noise (pial surface label) included in the ROIs from the total merged image (Array A1 below). This now should work fine as the calciumdxdettrial now correctly detects event onsets much more accurately than when this WaveCentroid script was written.

%James Ackman, 1/13/2011
%(1)Use similar dilation/erosion methods as for wave size calculation and make array of images masks for each wave
%(2)use image masks as input to regionprops to calculate weigthed centroids.

%Otsu's threshold for making wave center mask
% tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
% level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--

%Uncomment the plot commands at lines 155, 223, 228 to aid in troubleshooting.

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m . Use the with 'includebaseline = true' method.

if isfield(region,'waveonsets')
    
    %{
A = myMakeContourMovieWaves([],region,'true');
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
    
    %***Now setup contour image array A again from myMakeContourMovieWaves.m, this time with 'includebaseline = false' and pass this A array to the wavecentroid script below.
    A = myMakeContourMovieWaves([],region,'false');
    
    %--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
    c=1;
    f={};
    [szX,szY] = size(region.image);  %assuming szY is the largest dimension
    szZ = size(region.traces,2);
    if mod(max([szY szX]),min([szY szX])) == 0
        rXY=szY/szX;
        szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
    else
        rXY = 1;
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
    
    %-----New--this will get the histogram curve (trHist) for the movie that we will use to detect the half amplitude decrement for each wave
    % [hf1,trHist]=myPlotRasterHistManualPeaks(fnm,region,[],[],'true',spk,dec);
    % close(hf1)
    matAttOnOff=zeros(length(region.contours),length(region.traces));
    matAttOn=zeros(length(region.contours),length(region.traces));
    matAttOff=zeros(length(region.contours),length(region.traces));
    numCell=length(region.contours);
    for i=1:numCell
        ons=region.onsets{i};
        ofs=region.offsets{i};
        for j=1:length(ons)
            %         plot([ons(j) ofs(j)-1],[i i],'Linewidth',2)
            %         plot([ons(j)],[i ],'b.-')
            matAttOnOff(i,[ons(j):ofs(j)-1])=1;
            matAttOn(i,[ons(j)])=1;
            matAttOff(i,[ofs(j)])=1;
            %         covMat=xcov(region.traces(3,:)',sum(matAttOnOff)/numCell'
        end
    end
    trHist=sum(matAttOnOff)/numCell;
    %----------------
    
    region.wavecentroids={};
    wavecentroids={};
    for wfr=1:length(region.waveonsets)
%     for wfr=1:5
        disp([ 'wave no. ' num2str(wfr)])
        % wfr=8;
        
        halfampl = (trHist(region.wavepeaks(wfr)) - trHist(region.waveoffsets(wfr)))/2;
        halfindices=find(trHist(region.wavepeaks(wfr):region.waveoffsets(wfr)) < halfampl);
        halfamplindex = region.wavepeaks(wfr) + (halfindices(1) - 1);
        frames=region.waveonsets(wfr):halfamplindex;
        
        % frames=region.waveonsets(wfr):region.wavepeaks(wfr);
        % wavecentroid={}; %use structure if you want to keep more than one centroid for each waveframe
        wavecentroid=[]; %if you want to keep just the centroid of the max intensity wave area for each waveframe
        % BW_all=zeros(size(A(:,:,1)));
        BWmasks={};
        for d=1:length(frames)
            %         d=1;
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
            SE = strel('square', strel_sz(1)+3); %want close to a 3px radius to fill in the 1-2px gaps in between contours. e.g. 22x22 for 19x19 strel_sz ROIs
            % SE = strel('disk',5);
            %     A2=imclose(A1(:,:,fr),SE);  %this performs the sequential dilation and erosion operations
            % A2=imdilate(A2,SE);
            A2=imclose(A2,SE);
            %     figure;imshow(A2)
            
            level = graythresh(A2); %otsu's threshold method, including all zero pixels-- will separate into two populations
            % tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
            %     tmp=A1(:,:,wfr); tmp=tmp(:);  %Base otsu threshold on the total range of  maximum non-zero pixel values within the merged waveframe image
            %     level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--
            BW = im2bw(A2,level);
%                 figure; imshow(BW)
            [L, num] = bwlabel(BW, 8);
            % disp(num)
            RGB = label2rgb(L);
%             figure; imshow(RGB)
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
            %     STATS = regionprops(L,'Centroid');
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
                %         centroidlist=[centroidlist; STATS(i).Centroid];
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
            
            %     figure; imshow(BW)
            BW2=zeros(size(L));  %former
            [I,J]=find(L == idx(1));
            IND=sub2ind(size(L),I,J);
            BW2(IND)=1; %former
            %     BW_all(IND)=1;  %new
            BW2=im2bw(BW2);
            %     figure; imshow(BW2)
            
            %     BW3 = bwperim(BW2); %former
            %     %     BW3 = bwperim(BW_all); %new
            %     [row,col]=find(BW3);
            %     %     [I,J]=ind2sub(size(BW),find(BW3));
            %         figure; imshow(BW3)
            %     hold on;
            
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
            %{
    %the following is to limit the centroid frames to frames with non-overlapping centroids in the image masks
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
            %}
            
            disp(num2str(round(centroidlist)))
            
            %     disp(num2str(K))
            %here is final decision point
            %     if fr >= onset && fr<= frames(j)
            %     if fr >= onset && fr<= offset
            if K == 0
                %         figure; imshow(BW2)
                %         hold on;
                wavecentroid = [wavecentroid; centroidlist];
                %         plot(wavecentroid(d,1),wavecentroid(d,2),'+r')
                disp([num2str(fr) ' in, K = ' num2str(K)])
                BWmasks{length(BWmasks)+1}=BW2;
                %         hold off;
            else
                %         plot(centroidlist(1),centroidlist(2),'ob')
                wavecentroid=[wavecentroid; NaN NaN];
                disp([num2str(fr) 'out, K = ' num2str(K)])
                %         continue
            end
            % figure; plot(avgTrace)
            %----------------------------------------------------------------------
            
            % wavecentroid{fr}=[STATS.WeightedCentroid];
        end
        
        wavecentroids{wfr}=wavecentroid;  %save data in region structure array
        
        % %printout some results. If not enough detected waveframe centroids, no wavespeeds will be calculated.
        % ind=1:size(wavecentroid,1);
        % goodind=ind(find(~isnan(wavecentroid(:,1))));
        %
        % if length(goodind) > 1
        % diffind=diff(goodind)';
        % pxdist=sqrt((abs(diff(wavecentroid(goodind,1)))).^2 + (abs(diff(wavecentroid(goodind,2)))).^2);
        % % consecutivespeeds=(pxdist*region.spaceres)./(diffind*region.timeres)
        % % wavespeed=mean(consecutivespeeds)
        % % pxdist=pxdist(~isnan(pxdist));
        % else
        %     disp('not enough detected center of masses')
        % end
    end
    
    region.wavecentroids=wavecentroids;
    
else
    disp('run calciumdxDetectWaves.m first')
end