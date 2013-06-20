function [results] = makeBilateralWaveMaskMovie(wavesDBfnm,waveno1,waveno2)
%[results] = testXcorrDistMetric('wavesDB.mat',7,27);
%Peform 2D xcorr on pairs of waves isolated in waves database.
%James B. Ackman 2011-11-28

%waveno1 = 2; waveno2 = 11; wavesDBfnm = 'wavesDB.mat'; %TESTING
makefigs = 1;
makeCorr = 0;

fnmbase = ['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '_wv' num2str(waveno1) '-wv' num2str(waveno2)];  %base filename for saving graphics output
matObj = matfile(wavesDBfnm);  %use matfile to open connection to wavesDB without loading the whole thing into memory.
% szDB = size(matObj,'filename',2);

wave1data=matObj.imagedata(1,waveno1);   %load the data from the wavesDB. Matfile objects only uspport '()' indexing.
wave2data=matObj.imagedata(1,waveno2);

sz1 = size(wave1data.dfoverf); %we will set the no. of comparision frames to the length of the shorter wave array
sz2 = size(wave2data.dfoverf); %we will set the no. of comparision frames to the length of the shorter wave array
szZ = min([sz1(3) sz2(3)]);

A1 = wave1data.dfoverf(:,:,1:szZ);   %raw image data array for wave1
B1 = wave1data.wavemask(:,:,1:szZ);  %binary wave mask array for wave1
clear wave1data

%A2 = wave2data.dfoverf(:,:,1:szZ);   %raw image data array for wave2
B2 = wave2data.wavemask(:,:,1:szZ);  %binary wave mask array for wave2
clear wave2data

wave1orientation = matObj.orientation(1,waveno1);  %get the orientation values from DB for cropping
wave2orientation = matObj.orientation(1,waveno2);

wave1regionname = matObj.regionname(1,waveno1);  %get the region.names for determining which hemi images to flipud
wave2regionname = matObj.regionname(1,waveno2);

wave2regioncoords = matObj.regioncoords(1,waveno2);  %get the region.coords for determining size to crop template to

wave1maskparams = matObj.maskparams(1,waveno1);  %get the binary masking parameters
wave2maskparams = matObj.maskparams(1,waveno2);

%Test Xcorr2d between hemispheres
if makefigs > 0
    scrsz = get(0,'ScreenSize');
    fig1=figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(3)/2]);
    if makeCorr > 0
        fig2=figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(3)/2]);
    end
end
i=0;
dXY1 = [];
dXY2 = [];
maxCorr = [];
clear M1 M2
AlltextureDescriptors = [];
AllmomentInvariants = [];
Allstatxture = [];
Allregiondescriptors = [];

%	figure; imshow(regionMask1); 	figure; imshow(regionMask2);
cropSZ1=wave1orientation.orientation.value(1);  %so that we can crop the frame to just the hemisphere
cropSZ2=wave2orientation.orientation.value(1);

hbar2 = waitbar(0,'Running testXcorrDistMetric.m...');
for fr = 1:szZ
    waitbar(fr/szZ,hbar2);
    %fr = 87; %TESTING
    img1Raw = A1(:,:,fr);
    img1 = B1(:,:,fr);  %fetch the smoothed wavemask comparison frame
    img2 = B2(:,:,fr);
    %{
%---FOR raw dfoverf input with binary wavemask---
	img1 = A1(:,:,fr);  %fetch the raw data comparison frame
	img2 = A2(:,:,fr);
	
	img1 = gaussSmooth(img1,wave1maskparams.gaussSmoothSigma,'same');  %smooth the images by the same amount as was used for the getWaveMasks.m script
	img2 = gaussSmooth(img2,wave2maskparams.gaussSmoothSigma,'same');
    %}
    
    %prep coords for making an addition frame that is cropped to dimensions of region.coords to try as xcorr template--------------
    minX = min(wave2regioncoords.regioncoords(:,1));
    maxX = max(wave2regioncoords.regioncoords(:,1));
    minY = min(wave2regioncoords.regioncoords(:,2));
    maxY = max(wave2regioncoords.regioncoords(:,2));
    
    [m n] = find(img2);
    if ~isempty(n)
        minXtemplate = min(n);
        maxXtemplate = max(n);
    else
        minXtemplate = minX;
        maxXtemplate = maxX;
    end
    if ~isempty(m)
        minYtemplate = min(m);
        maxYtemplate = max(m);
    else
        minYtemplate = minY;
        maxYtemplate = maxY;
    end
    
    %------crop into two matching oriented images based on midline
    if strcmp(wave1regionname.regionname,'SC.R')  %assuming we want everything on the R-hemi coord system with anterior oriented towards the left of the image
        fr1=img1(1:cropSZ1,:);
        img1Raw=img1Raw(1:cropSZ1,:);
    elseif strcmp(wave1regionname.regionname,'SC.L')
        fr1=img1(cropSZ1+1:sz1(1),:);
        fr1=flipud(fr1);
        img1Raw=img1Raw(cropSZ1+1:sz1(1),:);
        img1Raw=flipud(img1Raw);
    else
        error('wrong regionnames, edit script to determine based on hemispheric coords...')
    end
    if strcmp(wave2regionname.regionname,'SC.R')  %assuming we want everything on the R-hemi coord system with anterior oriented towards the left of the image
        fr2=img2(1:cropSZ2,:);
        %		fr2crop = img2(minY:maxY,minX:maxX);
        fr2crop = img2(minYtemplate:maxYtemplate,minXtemplate:maxXtemplate);
    elseif strcmp(wave2regionname.regionname,'SC.L')
        fr2=img2(cropSZ2+1:sz2(1),:);
        %		fr2crop = img2(minY:maxY,minX:maxX);
        fr2crop = img2(minYtemplate:maxYtemplate,minXtemplate:maxXtemplate);
        fr2=flipud(fr2);
        fr2crop = flipud(fr2crop);
    else
        error('wrong regionnames, edit script to determine based on hemispheric coords...')
    end
    
    %	figure; imagesc(fr1)
    %	figure; imagesc(fr2)
    
    %{
%---FOR raw dfoverf input with binary wavemask----
	level = graythresh(img1);
	BW1 = im2bw(fr1,level);
	BW2 = im2bw(fr2,level);
	fr1ind = intersect(find(regionMask1crop),find(BW1));
	fr2ind = intersect(find(regionMask2crop),find(BW2));
	
	fr1bw = zeros(size(fr1));
	fr1bw(fr1ind) = fr1(fr1ind);  %set
	fr2bw = zeros(size(fr2));
	fr2bw(fr2ind) = fr2(fr2ind);
%	figure; imshow(fr1bw,[])
%	figure; imshow(fr2bw,[])
    %}
    fr1bw = fr1;
    fr2bw = fr2;
    
    [minDim,loca] = min([size(fr1,1) size(fr2,1)]);
    df = abs(diff([size(fr1,1) size(fr2,1)]));
    if loca == 1
        fr1bw = padarray(fr1bw, [df 0],'pre');
        %	img1Raw = padarray(img1Raw, [df 0],'pre');
    else
        fr2bw = padarray(fr2bw, [df 0],'pre');
    end
    
    img1Raw = padarray(img1Raw, [df 0],'pre');
    
    [minDim,loca] = min([size(fr1,2) size(fr2,2)]);
    df = abs(diff([size(fr1,2) size(fr2,2)]));
    if loca == 1
        fr1bw = padarray(fr1bw, [0 df],'post');
        %	img1Raw = padarray(img1Raw, [0 df],'post');
    else
        fr2bw = padarray(fr2bw, [0 df],'post');
    end
    
    img1Raw = padarray(img1Raw, [0 df],'post');
    
    if makefigs > 0
        disp(['padded fr1 size = ' num2str(size(fr1bw))]);
        disp(['padded fr2 size = ' num2str(size(fr2bw))]);
    end
    
    %	C = normxcorrn(fr1,fr2,'same');
    %	imagesc(C)
    %	M1(i) = getframe;
    %	[m n] = find(C == max(C(:)));
    %	dXY1 = [dXY1; m n];
    
    fr1bw = mat2gray(fr1bw);
    im1Raw = mat2gray(img1Raw);
    %	fr1bwNoise = imnoise(fr1bw,'gaussian');
    fr2bw = mat2gray(fr2bw);
    fr2bwNoise = imnoise(fr2bw,'gaussian');  %add gaussian background noise, so that corr value can always be calculated
    %	fr2bw = mat2gray(fr2bw);
    fr2crop = mat2gray(fr2crop); %for using the region.coords cropped frame as a template.
    
    %	if isempty(find(fr2crop))
    if isempty(find(fr1bw)) || isempty(find(fr2bw))
        %	fr2crop = imnoise(fr2crop,'gaussian');
        continue
    end
    
    
    i=i+1;  %frame counter for videoObject
    if makefigs > 0
        bothFrames = [fr1bw; ones(2,size(fr1bw,2)); flipud(fr2bwNoise)];
        set(0,'CurrentFigure',fig1)
        %	imagesc(bothFrames)
        %	colormap(gray)
        %	axis image
        
        imagesc(bothFrames)  %or img1Raw  %or fr1bw  %or bothFrames
        colormap(gray); axis image
        M1(i) = getframe;
    end
    template = fr2bwNoise;
    
    if makeCorr > 0
        %	C = normxcorr2(img1Raw,template);
        C = xcorrn(fr1bw,template,'same');
        %	C = xcorrn(fr1bw,fr2bwNoise,'same');
        %	C = xcorrn(fr1bw,fr2cropNoise,'same');
        if makefigs > 0
            set(0,'CurrentFigure',fig2)
            dfM = abs(diff([size(img1Raw,1) size(template,1)]));
            dfN = abs(diff([size(img1Raw,2) size(template,2)]));
            %	templateDisplay = padarray(template, [dfM dfN],1,'post');
            templateDisplay = padarray(C, [dfM dfN],0,'post');
            imagesc(templateDisplay); colormap(jet); %colorbar
            %	disp(num2str(size(C)))
            %	imagesc(templateDisplay)
            %	colormap(gray); axis image
            M2(i) = getframe;
        end
        
        mxVal = max(C(:));
        [m n] = find(C == mxVal);
        
        %	[mxVal, mxIdx] = max(abs(C(:)));
        %	[ypeak, xpeak] = ind2sub(size(C),mxIdx(1));
        %	corr_offset = [ (ypeak-size(template,1)) (xpeak-size(template,2)) ];
        %	m = mxIdx(1);
        
        if makefigs > 0
            disp(num2str(fr));
            disp(num2str(mxVal));
            disp(num2str(numel(m)));
        end
        if numel(m) > 1
            %	dXY2 = [dXY2; 1 1];
            continue
        else
            dXY2 = [dXY2; m n];
            %	dXY2 = [dXY2; corr_offset(1) corr_offset(2)];
        end
        maxCorr = [maxCorr; mxVal];
        
        
        %get some more image properties ast described on p. 651 of Gonzalez, Woods, Eddins Digital Image Processing with MATLAB
        %	textureDescriptors1 = getGraycoprops(fr1bw);
        texture1 = getStatxture(template);
        Allstatxture = [Allstatxture; texture1];
        
        textureDescriptors2 = getGraycoprops(template);
        AlltextureDescriptors = [AlltextureDescriptors; textureDescriptors2];
        
        momentInvariants = getInvmoments(template);
        AllmomentInvariants = [AllmomentInvariants; momentInvariants];
        
        %regionprops
        [L, num] = bwlabel(im2bw(template), 8);  %label connected regions
        STATS = regionprops(L, 'all'); %find areas
        areas=[STATS.Area]; %areas as single vector
        idx=find(areas == max(areas));
        idx = idx(1);
        
        regiondescriptors = [STATS(idx).Area STATS(idx).Centroid STATS(idx).Eccentricity STATS(idx).EquivDiameter STATS(idx).Extent STATS(idx).MajorAxisLength STATS(idx).MinorAxisLength STATS(idx).Orientation STATS(idx).Solidity];   %STATS(idx).EulerNumber leaving out because it gives a repeated value (no variation) and ends up causing a singularity in the pattern matrix when doing the matrix inverse calc for the Mahalanobis distance.
        Allregiondescriptors = [Allregiondescriptors; regiondescriptors];
    end
    %opticflow
end
close(hbar2)
if makefigs > 0
    close(fig1)
    if makeCorr > 0
        close(fig2)
    end
    %movie(M2,30)
end

%figure; mycolors = jet(size(dXY1,1));
%hold on;
%for i = 1:size(dXY1,1)
%plot(dXY1(i,2),dXY1(i,1),'o','Color',mycolors(i,:));
%end
%ylim([1 size(fr1,1)]); xlim([1 size(fr1,2)]); hold off;

if makeCorr > 0
    AllSqDistFromCenter = [];
    CenterXY = size(fr1bw) ./2;
end

if makefigs > 0
    if makeCorr > 0
        fig = figure;
        subplot(2,1,1)
        mycolors = jet(size(dXY2,1));
        hold on;
        
        for i = 1:size(dXY2,1)
            if makefigs > 0
                plot(dXY2(i,2),dXY2(i,1),'o','Color',mycolors(i,:));
                hold on
            end
            SqDistFromCenter = (dXY2(i,2) - CenterXY(1,2))^2 + (dXY2(i,1) - CenterXY(1,1))^2;
            AllSqDistFromCenter = [AllSqDistFromCenter; SqDistFromCenter];
        end
        if makefigs > 0
            ylim([1 size(fr1,1)]); xlim([1 size(fr1,2)]); hold off;
            axis square
            title(['xy xcorr displacements, ' 'wv' num2str(waveno1) '-wv' num2str(waveno2)])
            subplot(2,1,2)
            hist(maxCorr,20)
            axis square
            title('max xcorrn values')
            printfig('epsc',[fnmbase '.eps'],'false')
            close(fig)
        end
    end
end
%save([fnmbase '.mat'],'dXY2');

if makefigs > 0
    vidObj = VideoWriter([fnmbase '-wavemask' '.avi'])
    open(vidObj)
    for i =1:numel(M1)
        writeVideo(vidObj,M1(i))
    end
    close(vidObj)
    
    if makeCorr > 0
        vidObj = VideoWriter([fnmbase '.avi'])
        open(vidObj)
        for i =1:numel(M2)
            writeVideo(vidObj,M2(i))
        end
        close(vidObj)
    end
end

if makeCorr > 0
    results.dXY = dXY2;
    results.maxCorr = maxCorr;
    results.AllSqDistFromCenter = AllSqDistFromCenter;
    results.AlltextureDescriptors = AlltextureDescriptors;
    results.AllmomentInvariants = AllmomentInvariants;
    results.Allregiondescriptors = Allregiondescriptors;
    results.Allstatxture = Allstatxture;
end


function textureDescriptors = getGraycoprops(f)
%get some image texture properties for classification as described on p. 651 of Gonzalez, Woods, Eddins 'Digital Image Processing with MATLAB'
G2 = graycomatrix(f); %computes co-occurence matrix. Can also be called with more levels, G2 = graycomatrix(fr1bw,'NumLevels', 256)
G2n = G2/sum(G2(:)); %Normalized matrix
stats2 = graycoprops(G2,'all'); %Descriptors
maxProbability2 = max(G2n(:));
contrast2 = stats2.Contrast;
corr2 = stats2.Correlation;
energy2 = stats2.Energy;
hom2 = stats2.Homogeneity;
for I = 1:size(G2n,1);
    sumcols(I) = sum(-G2n(I,1:end).*log2(G2n(I,1:end) + eps));
end
entropy2 = sum(sumcols);
textureDescriptors = [maxProbability2 contrast2 corr2 energy2 hom2 entropy2];


function momentinvariants = getInvmoments(f)
%get some more image descriptors for classification as described on p. 660 of Gonzalez, Woods, Eddins 'Digital Image Processing with MATLAB'
phi = invmoments(f);
phinorm = -sign(phi).*(log10(abs(phi)));
momentinvariants = phinorm;


function texture1 = getStatxture(f)
%get some statistical image texture descriptors for classification based solely on the intensity histogram as described on p. 660 of Gonzalez, Woods, Eddins 'Digital Image Processing with MATLAB'
%returns row vector of 6 moments, mean, std dev, smoothness, third moment, uniformity, and entropy
t = statxture(f);
texture1 = t;