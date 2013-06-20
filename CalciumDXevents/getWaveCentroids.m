function region = getWaveCentroids(region)
%Used for calculating center of masses in sequential frames during a wave
%The resulting centroid list can be used in wave speed and wave direction calculations----

%revised--2/2/2011  Changed line 171 back to a gray threshold that simply excludes all non-zero pixels from the merged waveframe, instead of threshold based on intensites of all pixels in the image. This is because some images may have bright-non signal noise (pial surface label) included in the ROIs from the total merged image (Array A1 below). This now should work fine as the calciumdxdettrial now correctly detects event onsets much more accurately than when this WaveCentroid script was written.
%to include non-zero pixels, line 171 can be commented out and the previous or next lines can be used
%lines 292-301 can be commented out if overlapping centroid frames do not want to be commented out.
%updated to work by location, 5/6/2011

%James Ackman, 1/13/2011
%(1)Use similar dilation/erosion methods as for wave size calculation and make array of images masks for each wave
%(2)use image masks as input to regionprops to calculate weigthed centroids.

%Otsu's threshold for making wave center mask
% tmp=A2(:); %otsu's threshold method excluding zero pixels, to threshold wave center--
% level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--

%Uncomment the plot commands at lines 155, 223, 228 to aid in troubleshooting.

%***Must first setup up the contour image array A, from myMakeContourMovieWaves.m . Use the with 'includebaseline = true' method.

locationMarkers = unique(region.location);
regionsAll = splitRegion(region);
for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    wavecentroids = getWaveCentroidsByLocation(tmpregion,locationIndex);
    region.wavedata{locationIndex}.wavecentroids=wavecentroids;
end

function wavecentroids = getWaveCentroidsByLocation(region,locationIndex)
if isfield(region.wavedata{locationIndex},'waveonsets')
    
    % A = myMakeContourMovieWaves([],region,'true');
    %{
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
    figure; imshow(Adummy)
    %Add dummy frame to master wave list length of region.wavedata{locationIndex}.waveonsets
    A1(:,:,i)=Adummy;
end
    %}
    
    %***Now setup contour image array A again from myMakeContourMovieWaves.m, this time with 'includebaseline = false' and pass this A array to the wavecentroid script below.
    A = myMakeContourMovieWavesOnsets([],region,locationIndex,'false');
%     A = myMakeContourMovieWaves([],region,'false');
    strel_sz = [];
    for i = 1:length(region.contours)
        out = getStrelsize(region,i);
        strel_sz = [strel_sz; out];
        disp(i)
    end
    strel_sz=round(mean(strel_sz,1));
    
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
    
    %trace Histogram with waveonsets and offsets for testing. Same plot from calciumdxDetectWaves.m---------------
    %{
    figure;
    a(1) = subplot(2,1,1);
    plot(trHist);
    hold on;
    plot(region.wavedata{locationIndex}.wavepeaks,trHist(region.wavedata{locationIndex}.wavepeaks),'or');
    plot(region.wavedata{locationIndex}.waveonsets,trHist(region.wavedata{locationIndex}.waveonsets),'ok');
    plot(region.wavedata{locationIndex}.waveoffsets,trHist(region.wavedata{locationIndex}.waveoffsets),'og');
    hold off;
    a(2) = subplot(2,1,2);
    imagesc(dfoverf(region.traces));
    linkaxes(a,'x')
    %}
    
    %----------------------------------------------------------------------
    %use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19px contours.
%     SE = strel('rectangle', [2 2]);
%     A2=imdilate(A1(:,:,wfr),SE);
%     SE = strel('square', strel_sz(1)+3); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
%     A2=imclose(A2,SE);
%     
%     SE = strel('disk',10);
% %     A2=imdilate(A1(:,:,wfr),SE);
%     A2=imdilate(A2,SE);
%     SE = strel('disk',2);
% %     A2=imdilate(A1(:,:,wfr),SE);
%     A2=imclose(A2,SE);
%     
%     figure;imshow(A2)
% %     figure; contour(interp2(A2,1));
%     
%     [X,Y]=meshgrid(1:size(A2,2),1:size(A2,1));
%     [DX,DY] = gradient(A2,5,5);
%     figure; contour(X,Y,A2)
%     hold on
%     quiver(X,Y,DX,DY)
% %     colormap hsv
%     hold off
%     figure; imagesc(DX); figure; imagesc(DY)
    
    
    %----------------------------------------------------------------------
    hbar = waitbar(0,'Please wait...');
    %     region.wavedata{locationIndex}.wavecentroids={};
    wavecentroids={};
    for wfr=1:length(region.wavedata{locationIndex}.waveonsets)
        %     for wfr=1:5
        disp([ 'wave no. ' num2str(wfr)])  %for testing
        % wfr=8;
        
        %         %plot gray dotted outline of each region----
        %         figure;
        %         imshow(zeros(size(region.image)));
        %         hold on
        %         for numcoords = 1:length(region.coords)
        %             plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
        %         end
        %         hold off
        
        %find the 50% decay time-------------------------------------------
%         halfampl = (trHist(region.wavedata{locationIndex}.wavepeaks(wfr)) - trHist(region.wavedata{locationIndex}.waveoffsets(wfr)))/2;
%         %         disp(num2str(halfampl))
%         halfindices=find(trHist(region.wavedata{locationIndex}.wavepeaks(wfr):region.wavedata{locationIndex}.waveoffsets(wfr)) < halfampl+trHist(region.wavedata{locationIndex}.waveoffsets(wfr)));
%         %                 disp(num2str(halfindices))
%         halfamplindex = region.wavedata{locationIndex}.wavepeaks(wfr) + (halfindices(1) - 1);
        %                 disp(num2str(halfamplindex))
%         frames=region.wavedata{locationIndex}.waveonsets(wfr):halfamplindex;
%                 frames=region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.wavepeaks(wfr);
                frames=region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.waveoffsets(wfr);
%         disp(num2str(frames))
        
        % frames=region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.wavepeaks(wfr);
        % wavecentroid={}; %use structure if you want to keep more than one centroid for each waveframe
        wavecentroid=[]; %if you want to keep just the centroid of the max intensity wave area for each waveframe
        % BW_all=zeros(size(A(:,:,1)));
%         BWmasks={};

        szA = size(A(:,:,1));
        BWmasks = zeros([szA length(frames)]);
        for d=1:length(frames)
%             d=48;
            fr=frames(d);
            %     wavecentroid{d}=[];
%                 fr=44; %testing with diff frames
            % figure;imshow(A1(:,:,fr))
            %want close to a 3px radius to fill in the 1-2px gaps in between contours.
            SE = strel('rectangle', [2 2]);   %<--------------------
            % SE = strel('disk',5);
            % A2=imdilate(A1(:,:,wfr),SE);
            A2=imdilate(A(:,:,fr),SE);  %<----------------------
%                         figure;imshow(A2)  %TESTING
%             figure; imshow(A(:,:,fr))  %TESTING
            
            %try imclose
            %use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  e.g. 22x22 for 19x19px contours. Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19 strel_sz ROIs
            SE = strel('square', strel_sz(1)+3); %<----------------------
            %     A2=imclose(A1(:,:,fr),SE);  %this performs the sequential dilation and erosion operations     
            A2=imclose(A2,SE);  %<----------------------
            BW = im2bw(A2);  %make binary
            [L, num] = bwlabel(BW, 8);  %label connected regions
            RGB = label2rgb(L); %for display testing
%             figure; imshow(RGB); %for display TESTING
            STATS = regionprops(L, 'Area'); %find areas
            areas=[STATS.Area]; %areas as single vector
            minIdx = find(areas <= 1*prod(strel_sz+1));  %get labeled region indices that are less than or equal to one strel_sz area with the 1 pixel border added by imdilate above
            if ~isempty(minIdx)
                for minIdx1 = minIdx
                   A2(L==minIdx1) = 0;
                end
            end
            
            %----dilate connected regions to help with centroid calculations---
%             SE = strel('disk',10);
            SE = strel('disk',round(strel_sz(1)*2.5));  %2.5x bigger than strel_sz. i.e. for 4px square strel that give approx 20um diam regions, this will give 10px diam that will give a 50-80um region.
            A2=imdilate(A2,SE);
            SE = strel('disk',2);
            A2=imclose(A2,SE);
%             figure;imshow(A2)  %TESTING

%             SE = strel('square', 4*strel_sz(1));
%             A2 = imerode(A2,SE);

            level = graythresh(A2); %otsu's threshold method, including all zero pixels-- this will simply separate into two populations
%             tmp=A2(:); %otsu's threshold method excluding zero pixels , to threshold wave center-- Good for CCD recordings and such that have high pixel densities and lower noise.
%                             tmp=A1(:,:,wfr); tmp=tmp(:);  %Base otsu threshold on the total range of  maximum non-zero pixel values within the merged waveframe image.  Can be useful for noisier movies with less pixel density or ones where the signal detection procedure is in question. Not useful for CCD based recordings.
%             level=graythresh(tmp(tmp>0));  %otsu's threshold method excluding zero pixels, to threshold wave center--
            BW = im2bw(A2,level);
%             figure; imshow(BW)   %TESTING
            [L, num] = bwlabel(BW, 8);
            % disp(num)
            RGB = label2rgb(L);
%             figure; imshow(RGB)  %TESTING
            
            if num == 0
                %         wavecentroid{d}=[wavecentroid{d}; NaN];
                wavecentroid=[wavecentroid; NaN NaN];
%                 disp([num2str(fr) 'out first'])  %for testing
                continue %continues onto next iteration of for loop
            end
            
            %=======Limit regions for centroid analysis by max area=======================================
            %{
            STATS = regionprops(L, 'Area');
            areas=[STATS.Area];
            idx=find(areas == max(areas));
            
            % disp(fr);
            % disp([STATS.Area])
            
            if length(idx) > 1  %in case there are multiple detected areas of same size, like 3x size(ROI), then probably means the detected transients in the wave are too noisy to reliably detect the waveframe center anyways.
                wavecentroid=[wavecentroid; NaN NaN];
                disp([num2str(fr) 'out next'])  %for testing
                continue
            end
            
            if areas(idx(1)) <= 3*prod(strel_sz+1)  %strel_sz is the the dimensions -1px of the rois used in region.contours for analysis. This gets rid of ones that are not big enough to be part of a wave for detecting the true wave center of mass.
                %         wavecentroid{d}=[wavecentroid{d}; NaN];
                wavecentroid=[wavecentroid; NaN NaN];
                disp([num2str(fr) 'out second'])  %for testing
                continue %continues onto next iteration of for loop
            end
            %}
            %========Perform centroid calculation based on labeled regions================================
            %         figure; imshow(RGB)
            %     hold
%             STATS = regionprops(L,A2,'WeightedCentroid');
                STATS = regionprops(L,'Centroid');
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
            idx = 1:length([STATS.Centroid])/2;
            centroidlist=[];
            for i=idx
%                 centroidlist=[centroidlist; STATS(i).WeightedCentroid];
                        centroidlist=[centroidlist; STATS(i).Centroid];
            end
            meanCentroid = mean(centroidlist,1);
%             medianCentroid = median(centroidlist,1);
%             figure; imshow(A2); hold on;
%             plot(centroidlist(:,1),centroidlist(:,2),'.b')
%             plot(meanCentroid(1),meanCentroid(2),'xg')
%             plot(medianCentroid(1),medianCentroid(2),'xr')
%             hold off;
            
            %----------------------------------------------------------------------
            %algorithm for limiting the moving centroid frames to those in the wavefront
            % CC = bwconncomp(BW,8);
            % numPixels = cellfun(@numel,CC.PixelIdxList);
            % [biggest,idx] = max(numPixels);
            % BW2=zeros(size(BW));
            % BW2(CC.PixelIdxList{idx}) = 1;
            % figure; imshow(BW)
            % figure; imshow(BW2)
            
            %======Get binary image of max area region used for centroid calc=========================
            %-----used in the code below for setting the K parameter for determiing if the centroid frame had been obtained yet or whether it should be skipped. 
            %     figure; imshow(BW)
            BW2=zeros(size(L));  %former
            [I,J]=find(L == idx(1));   %get indices for biggest labeled region in L from regionprops
            IND=sub2ind(size(L),I,J);
            BW2(IND)=1; %former  %set BW2 mask to the mask for biggest labeled region in L  from regionprops
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
        id=find(region.onsets{inp(c)} >= region.wavedata{locationIndex}.waveonsets(wfr) & region.onsets{inp(c)} <= region.wavedata{locationIndex}.waveoffsets(wfr));
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
            %=====================================================================
            %the following is to limit the centroid frames to frames with non-overlapping centroids in the image masks
            %{
            if ~isempty(BWmasks)
                for i = 1:length(BWmasks)
                    %             inp=inpolygon(meanCentroid(1),meanCentroid(2),BWmasks{i}(:,2),BWmasks{i}(:,1));  %find those that are bounded by the original ROI image mask
                    if BWmasks{i}(round(meanCentroid(2)),round(meanCentroid(1)))
                        K=1;
                    end
                    %             if inp == 1
                    %                 %     disp(num2str(K))
                    %                 K=1;
                    %             end
                end
            end
            %}
            %=====================================================================
            %             disp(num2str(round(meanCentroid)))  %for testing
            
            %     disp(num2str(K))
            %here is final decision point
            %     if fr >= onset && fr<= frames(j)
            %     if fr >= onset && fr<= offset
            if K == 0
                wavecentroid = [wavecentroid; meanCentroid];
%                 BWmasks{length(BWmasks)+1}=BW2;
                %-----------for TESTING------------------------------------
                %{
%                 figure; imshow(RGB)
                figure; imshow(BW2)
%                 figure;imshow(A2)
                hold on
                for numcoords = 1:length(region.coords)
                    plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
                    plot(wavecentroid(d,1),wavecentroid(d,2),'+r')  %for TESTING
                    %                 disp([num2str(fr) ' in, K = ' num2str(K)]) %for testing
                end
                hold off
                %-----------END for TESTING--------------------------------
                %}
            else
                wavecentroid=[wavecentroid; NaN NaN];
                %                 hold on; %for TESTING
                %                         plot(meanCentroid(1),meanCentroid(2),'ob') %for TESTING
                %                 hold off; %for TESTING
                %                 disp([num2str(fr) 'out, K = ' num2str(K)])  %for testing
                %         continue
                
            end
            % figure; plot(avgTrace)
            %----------------------------------------------------------------------
            BWmasks(:,:,d) = A2;
            % wavecentroid{fr}=[STATS.WeightedCentroid];
        end
        
        disp(['centroid frames: '])
        disp(num2str(wavecentroid))
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
        waitbar(wfr/length(region.wavedata{locationIndex}.waveonsets),hbar);
    end
    
    %     region.wavedata{locationIndex}.wavecentroids=wavecentroids;
    
else
    disp('run calciumdxDetectWaves.m first')
end
close(hbar)
