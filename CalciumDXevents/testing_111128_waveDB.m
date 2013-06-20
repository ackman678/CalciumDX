wavesDBfnm = 'wavesDB_v2-raw.mat';
wavesDBfnm = 'wavesDB_v2-smooth.mat';
wavesDBfnm = 'wavesDB.mat';
wavesDBfnm = 'wavesDB_tmp.mat';

matObj = matfile(wavesDBfnm);
szDB = size(matObj,'filename',2);
allOns = [region.wavedata{2}.waveonsets region.wavedata{3}.waveonsets];
allOffs = [region.wavedata{2}.waveoffsets region.wavedata{3}.waveoffsets];

for i = 1:szDB
    tmp=matObj.imagedata(1,i);
    sz1 = size(tmp.dfoverf);
    sz2 = size(tmp.wavemask);
    disp(['wave no. ' num2str(i)])
    disp(['sz1: ' num2str(sz1)])
    disp(['sz2: ' num2str(sz2)])
    disp(['wvSz: ' num2str(allOffs(i) - allOns(i)) ' frames'])
end


%---take a look at the array masks for a fetched wave. -----
tmp=matObj.imagedata(1,423); 
implay(mat2gray(tmp.dfoverf)); 
implay(mat2gray(tmp.wavemask));




%-----Batch testXcorrMetric on multiple waves in a hemisphere------
clear data;
wave1number = 121;
% wave2numbers = [1:18];
% wave2numbers = [1:59];
wave2numbers = [14 121];
% wave2numbers = [47:59];
waveCorrRefno = 14;
allmeanDist = [];
allmedianDist = [];
h = waitbar(0,'Please wait...');
for i=wave2numbers
    waitbar(find(wave2numbers==i)/numel(wave2numbers),h);
    results = testXcorrDistMetric('wavesDB.mat',wave1number,i);
    data(i).xcorrResults = results;
    
    %{
    disp(['wave' num2str(wave1number) ' - wave' num2str(i) ' mean distance'])
    %     meanDist = round(mean(data(i).AllSqDistFromCenter));
    
    dXY2 = data(i).dXY;
    AllSqDistFromCenter = [];
    for j = 1:size(data(i).dXY,1)
        CenterXY = data(waveCorrRefno).dXY(j,:);
        SqDistFromCenter = (dXY2(i,2) - CenterXY(1,2))^2 + (dXY2(i,1) - CenterXY(1,1))^2 + (data(i).maxCorr(j) - data(waveCorrRefno).maxCorr(j))^2;
        AllSqDistFromCenter = [AllSqDistFromCenter; SqDistFromCenter];
    end
    meanDist = round(mean(AllSqDistFromCenter));
    
    allmeanDist = [allmeanDist; meanDist];
    disp(num2str(meanDist))
    medianDist = round(median(data(i).AllSqDistFromCenter));
    allmedianDist = [allmedianDist; medianDist];
    disp(num2str(medianDist))
    
    %}
end
close(h)
figure; plot(wave2numbers,allmeanDist,'-ok'); title(['meanDist, wave' num2str(wave1number) ' vs random'])
figure; plot(wave2numbers,allmedianDist,'-ok'); title(['medianbDist, wave' num2str(wave1number) ' vs random'])

% %get max dX and max dY values to normalize dXY and test the scale invariance of mahal dist calc
% x = [];
% y = [];
% for i = wave2numbers
%     x = [x; max(data(i).dXY(:,1))];
%     y = [y; max(data(i).dXY(:,2))];
% end
% maxX = max(x);
% maxY = max(y);


allmeanDistMahal  = [];
allmeanEuclDist = [];
nRows = size(data(waveCorrRefno).dXY,1);  %for using observed reference wave as reference
randJitter = rand(nRows,2);  %for using observed reference wave as reference
dXY1 = data(waveCorrRefno).dXY + randJitter;  %for using observed reference wave as reference
for i = wave2numbers
    disp(num2str(i))
%     i=27
%     X = [data(waveCorrRefno).dXY data(waveCorrRefno).maxCorr]; %reference sample
%     Y = [data(i).dXY data(i).maxCorr];  %observation sample

%     X = [data(waveCorrRefno).dXY data(waveCorrRefno).AlltextureDescriptors data(waveCorrRefno).AllmomentInvariants]; %reference sample
%     Y = [data(i).dXY data(i).AlltextureDescriptors data(i).AllmomentInvariants];  %observation sample
%     X = [data(waveCorrRefno).Allregiondescriptors data(waveCorrRefno).dXY data(waveCorrRefno).maxCorr data(waveCorrRefno).AlltextureDescriptors data(waveCorrRefno).AllmomentInvariants data(waveCorrRefno).Allstatxture]; %reference sample
%     Y = [data(i).Allregiondescriptors data(i).dXY data(i).maxCorr data(i).AlltextureDescriptors data(i).AllmomentInvariants data(i).Allstatxture];  %observation sample
    
    regpropInd = [1:10];  %in case you want to limit the no. of properties input to the calculation from regionprops
    textureInd = [1:2 4:6];
%     X = [data(waveCorrRefno).dXY data(waveCorrRefno).maxCorr data(waveCorrRefno).AlltextureDescriptors data(waveCorrRefno).AllmomentInvariants data(waveCorrRefno).Allstatxture data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
%     Y = [data(i).dXY data(i).maxCorr data(i).AlltextureDescriptors data(i).AllmomentInvariants data(i).Allstatxture data(i).Allregiondescriptors(:,regpropInd)];  %observation sample    
   
%     dXY1 = [data(waveCorrRefno).dXY(:,1)./maxX data(waveCorrRefno).dXY(:,2)./maxY];  %just to test the scale invariance of mahal
%     dXY2 = [data(i).dXY(:,1)./maxX data(i).dXY(:,2)./maxY];

    if waveCorrRefno == wave1number
        X = [dXY1 data(waveCorrRefno).AlltextureDescriptors(:,textureInd) data(waveCorrRefno).AllmomentInvariants data(waveCorrRefno).Allstatxture data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
    else
        X = [data(waveCorrRefno).dXY data(waveCorrRefno).AlltextureDescriptors(:,textureInd) data(waveCorrRefno).AllmomentInvariants data(waveCorrRefno).Allstatxture data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
%         X = [data(waveCorrRefno).dXY data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
    end

    if ~isempty(data(i).dXY)
    Y = [data(i).dXY data(i).AlltextureDescriptors(:,textureInd) data(i).AllmomentInvariants data(i).Allstatxture data(i).Allregiondescriptors(:,regpropInd)];  %observation sample            
%     Y = [data(i).dXY data(i).Allregiondescriptors(:,regpropInd)];  %observation sample            
    %also  data(i).maxCorr  data(i).AllSqDistFromCenter
%     RefDistMahal = mahal(X,X);   %for ttest2
    DistMahal = mahal(Y,X); % Mahalanobis
    else
        DistMahal = 1e05;  %arbitrary large Mahal distanc number. e.g. if there were no frame information returned from the current comparision wave 
    end
    

    dXY2 = data(i).dXY;
    AllSqDistFromCenter = [];
    for j = 1:size(dXY2)
        CenterXY = data(waveCorrRefno).dXY(1,:);
        SqDistFromCenter = (dXY2(j,2) - CenterXY(1,2))^2 + (dXY2(j,1) - CenterXY(1,1))^2;
        AllSqDistFromCenter = [AllSqDistFromCenter; SqDistFromCenter];
    end
    meanDist = round(mean(AllSqDistFromCenter));
    
    allmeanEuclDist = [allmeanEuclDist; meanDist];
    
%         [h,p] = ttest2(RefDistMahal,DistMahal)  %not useful, everything is significant
%         DistMahal = mahalanobis(Y,X); % Mahalanobis
%     DistMahal = pdist2(Y,X,'euclidean'); % Mahalanobis
%     D = sqrt(sum(abs(bsxfun(@minus,permute(X,[1 3 2]), permute(Y,[3 1 2]))).^2,3)); %Euclidean distance between pattern vector populations
%     meanEuclDist = mean(D);
%     allmeanEuclDist = [allmeanEuclDist; meanEuclDist];

%         figure; hist(DistMahal,20)
%         figure; plot(1:length(DistMahal),DistMahal,'-ok')
%     figure; imagesc(DistMahal); colorbar; title([num2str(waveCorrRefno) ' - ' num2str(i)])
    meanDistMahal = nanmean(DistMahal);
%     disp(num2str(meanDistMahal));
    allmeanDistMahal = [allmeanDistMahal; meanDistMahal];
%     allmedDistMahal = [allmedDistMahal; medDistMahal];
end
figure; plot(wave2numbers,allmeanDistMahal,'-ok'); title(['mean Mahalanobis dist from ref wave' num2str(waveCorrRefno) ' vs random']); zoom yon

figure; plot(wave2numbers,allmeanEuclDist,'-ok'); title(['Eucl dist from ref wave' num2str(waveCorrRefno) ' vs random']); zoom yon; 

save('xcorrn_testing_111130_wv6.mat','dataWave6')
save('xcorrn_testing_111130_wv14-34_110323_08.mat','dataWave14')

getWaveMasks({filename},region,[],series1);  %for a single file

% X = mvnrnd([0;0],[1 .9;.9 1],100); %mahal demo
% Y = [1 1;1 -1;-1 1;-1 -1]; %mahal demo
% d1 = mahal(Y,X)


[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(data,wave1number,wave1number,waveCorrRefno);
[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(data,wave1number,wave2numbers(1),waveCorrRefno);

/Users/ackman/Data/2photon/110308i/110308_05_defaultROIs_kalman2.mat

%-----Get list for indices for wavesDB------
wavesDBfnm = 'wavesDB.mat';
% wavesDBfnm = 'wavesDB_tmp.mat';
matObj = matfile(wavesDBfnm);
szDB = size(matObj,'filename',2);
% for i = 1:szDB
%     fnm=matObj.filename(1,i);
%     wvnum=matObj.wavenumber(1,i);
%     disp([num2str(i) ' ' fnm.filename ' ' num2str(wvnum.wavenumber)])
% end

%the following is a much faster, more efficient way of getting the wavesDB lookup table.
waveDBidx = [1:szDB]';
tmp = matObj.filename;
filename = {tmp.filename}';
tmp = matObj.regionname;
regionname = {tmp.regionname}';
tmp = matObj.wavenumber;
wavenumber = {tmp.wavenumber}';
for i = 1:szDB
    disp([num2str(waveDBidx(i)) ' ' filename{i} ' ' num2str(regionname{i}) ' ' num2str(wavenumber{i})])
end





%-------Reopen image and test new avgimg for dfoverf array, and median filtering for making image mask---
[data, series1] = myOpenOMEtiff([]);
avgIm = mean(series1(:,:,448:468),3);
waveONidx = 485; waveOFFidx = 585;
A = double(series1(:,:,waveONidx:waveOFFidx));  %make raw image array
sz = size(A);  %array size

for fr = 1:sz(3)
    A(:,:,fr) = (A(:,:,fr) - avgIm)./avgIm;   %replace each raw image with the normalized dF/F representation
end

BWarray1 = [];
gaussSmoothSigma = 3;
graythreshlevel = [];
absoluteMinimumDfoverfValue = [];
for fr = 1:sz(3)
    img1 = gaussSmooth(A(:,:,fr),gaussSmoothSigma,'same');  %accentuate big blobby objects at 6 sigma, about 50-100um diam.
    
    %               fr1=img1(1:cropSZ,:);   %use these lines for cropping and flipping the masks for comparisons between hemispheres.
    %                fr2=img1(cropSZ+1:sz(1),:);
    %                fr2=flipud(fr2);
    
    level = graythresh(img1);
    BW = im2bw(img1,level);
    BWarray1(:,:,fr) = BW;
end
implay(BWarray1)
tmp=matObj.imagedata(1,2); implay(mat2gray(tmp.dfoverf)); implay(mat2gray(tmp.wavemask));







%test normxcorr2 example
%   Example
%   -------
template = .2*ones(11); % make light gray plus on dark gray background
template(6,3:9) = .6;
template(3:9,6) = .6;
BW = template > 0.5;      % make white plus on black background
figure, imshow(BW), figure, imshow(template)
% make new image that offsets the template
offsetTemplate = .2*ones(101);
offset = [70 70];  % shift by 3 rows, 5 columns
offsetTemplate( (1:size(template,1))+offset(1),...
    (1:size(template,2))+offset(2) ) = template;
figure, imshow(offsetTemplate)
%
% cross-correlate BW and offsetTemplate to recover offset
cc = normxcorr2(BW,offsetTemplate);
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
corr_offset = [ (ypeak-size(template,1)) (xpeak-size(template,2)) ];
isequal(corr_offset,offset) % 1 means offset was recovered

