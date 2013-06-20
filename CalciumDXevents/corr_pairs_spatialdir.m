function results = corr_pairs_spatialdir(region,numres,fnm,locationIndexPair)
%corr_pairs_spatialdir.m
%results = corr_pairs_spatialdir(region,1000,fnm)
%identify correlated ROIs which have correlated spike times, spatial locations, and wave directions
%created by James B. Ackman 2011-04-22
close('all')
if nargin < 4 || isempty(locationIndexPair), locationIndexPair = [2 3]; end

% load matlab3.mat
% load corrspatialtest_111108.mat

fname2base = [fnm(1:end-4) 'dCorrTDL_'];


% numCorrThreshold = 2;  %No. of temporal-correlations that a cell pair must reach to be considered a significant correlation for counting obs and res TDL corrs. Define for specific hypothesis testing or define automatically in relation to observed no. of spatio-temporal-directive correlations found below at script starting on line 97
sig = 1;
nDir = 4;
% roiBin_width = 25  %set automatically to 50 or 25px for 5x or 2.5x CCD movies in makeMatchingROImeshgrid.m

%?CHANGE-- need setup variables of no. of possible locations and directions and change any instance of 20 or 4 below.
% sz = size(region.traces,2);
% s = zeros(length(region.onsets),sz);
% for i = 1:size(region.traces,1)
%     s(i,region.onsets{i}) = 1;
% end


%---Create spike array movie---
A = myMakeContourMovieWavesOnsets([],region);

%---Dilate ROIs to connect neighbors and exclude isolated pixels having no neighbors nearby (activations not within the activity wave front)
strel_sz = getStrelsize(region);
[A2,A3] = dilateContours(A,strel_sz);

%---Setup new spatially binned ROI spike array containing spikes from array found from the dilated spike array A2
[bigmesh] = makeMatchingROImeshgrid(region,locationIndexPair);
s = intersectSpikeArrayContours(bigmesh.contours,A2);
%}

%{
%set spike array to just first spike for each cell inside each wave.--------------
ntSpk = zeros(size(s));
for numCell = 1:size(s,1)
    locationIndex = bigmesh.location(numCell);
    spks = find(s(numCell,:));
    if ~isempty(spks)
    for d = 1:numel(region.wavedata{locationIndex}.waveonsets);
        ind = find(spks >= region.wavedata{locationIndex}.waveonsets(d) & spks <= region.wavedata{locationIndex}.waveoffsets(d));
        if ~isempty(ind)
        ntSpk(numCell,spks(ind(1))) = 1;
        end
    end
    end
end
s = ntSpk;
%}

for c = 1:size(s,1)
    for findSpks = find(s(c,:))
        in = find(abs((1:size(s,2))-findSpks)<=1);  %smooth the spikes with +/- 1 signal
        s(c,in) = s(c,findSpks);
    end
end
figure; imagesc(s); title('s'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
%--------------------------------------------------------------------------------

%--Setup location spike matrix----------------
%L = region.locations  %need to setup
L = setupLocationMatrix(s,bigmesh,region);

%--Setup direction spike matrix----------------
%D = region.directions  %need to setup
D = setupDirectionMatrix(s,bigmesh,region);

figure; imagesc(L); title('L'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(D); title('D'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
%}

%{
%----RANDOM spike matrices for TESTING--------------------------
% results = [];
% for ct = 1:100
test.numCells = 2;
test.numCorr = 2;
test.numSpikes = 10;
test.numCorrSpikes = 2;
s = zeros(test.numCells,size(region.traces,2));
for i=1:size(s,1)
    c1 = monte(test.numSpikes,size(s,2));
    s(i,c1) = 1;
end
test.CorrIdx = fix(rand(1,test.numCorrSpikes)*size(s,2))+1;
s(1:test.numCorr,test.CorrIdx) = 1;

L = zeros(test.numCells,size(region.traces,2));
D = zeros(test.numCells,size(region.traces,2));
L(s>0) = fix(rand(1,numel(find(s>0)))*20)+1;  %20 different ROI simulated locations
D(s>0) = fix(rand(1,numel(find(s>0)))*4)+1;  %4 different cardinal directions

test.CorrLoca = fix(rand(1,length(test.CorrIdx))*20)+1;
test.CorrDir = fix(rand(1,length(test.CorrIdx))*4)+1;

for i=1:length(test.CorrLoca)
    L(1:test.numCorr,test.CorrIdx(i)) = test.CorrLoca(i);
    D(1:test.numCorr,test.CorrIdx(i)) = test.CorrDir(i);
end
%}

%{
%setup IMAGES of spike arrays for printing and SMOOTH s,L,D ---------------------------------
img.s = zeros(size(s));
img.L = zeros(size(s));
img.D = zeros(size(s));
for c = 1:size(s,1)
    for f = find(s(c,:))
        in = find(abs((1:size(s,2))-f)<=1);  %smooths the spikes +/-1 frame for pretty printing
        img.s(c,in) = s(c,f);
        img.L(c,in) = L(c,f);
        img.D(c,in) = D(c,f);
    end
end
figure; imagesc(img.s); title('s')
figure; imagesc(img.L); title('L')
figure; imagesc(img.D); title('D')
%----END example spike matrices for testing---------------------
%}

%---Setup img spike arrays without extra smoothing when using s array based on output from myMakeContourMoviesWavesOnsets.m
img.s = s;
img.L = L;
img.D = D;

%---Get observed correlation values---
% s_thick = gauss_events(s,sig);
s_thick = s; %if the spike matrix has effectively already been smoothed spatially and temporally from myMakeContourMoviesWavesOnsets.m
corr = s_thick*s_thick';
% corr = s*s';  %only if you will test without gaussian spike smoothing

%---TESTING to count num observed corr spikes. **uses img.L, and img.D above**--------------
[r,p] = corrcoef(s_thick');
tmp = tril(p,-1);
tmp = tmp + triu(ones(size(tmp)));
nanvalues = isnan(tmp);
tmp(nanvalues) = 1;

figure; imagesc(tmp); title('tmp'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])

[iPairs jPairs] = find(tmp < 0.05);
obs.countn = size(iPairs,1);
obs.countL = zeros(size(corr));
obs.countD = zeros(size(corr));
obs.count = zeros(size(corr));  %setup the results histogram array
numCorrThreshold=1;

    if ~isempty(iPairs)
        for c=1:size(iPairs,1)
            %get the corr cell indices i,j
            %         size([i j],1);
            i = iPairs(c);
            j = jPairs(c);
            %find the correlated spike indices for cell i, cell j
            cellIND = s_thick(i,:) .* s_thick(j,:);  %element wise multiplication to intersect the spike times for the two corr cells
            %             level = graythresh(cellIND);  %separate the histrogram of corr values into two populations, especially useful if gauss events was used for spike smoothing and there is loads of near-zero values.
            [muhat, sigmahat] = normfit(cellIND);
            level = muhat+(2*sigmahat);
            spkIND = find(cellIND>level);
            dfSpk = find(diff(spkIND) > 1);
            if ~isempty(dfSpk)
%                 spkIND = [spkIND(dfSpk) spkIND(dfSpk+1)];
                spkIND = setdiff(spkIND,spkIND+1);
            end
            %--------------------------------------------------------------
            %test spatial match and direction match at the correlated spike indices in smoothed L and D
            ind.L.i = img.L(i,spkIND);
            ind.L.j = img.L(j,spkIND);
            ind.D.i = img.D(i,spkIND);
            ind.D.j = img.D(j,spkIND);
            %numel n. overlapping spikes
            intL = ismember(ind.L.i(ind.L.i>0),ind.L.j(ind.L.j>0)); %fyi, intersect sorts the resulting set in ascending order.  20110901 changed to ismember from intersect, as intersect will not return repeated values
            intD = ismember(ind.D.i(ind.D.i>0),ind.D.j(ind.D.j>0)); %fyi, intersect sorts the resulting set in ascending order   20110901 changed to ismember from intersect, as intersect will not return repeated values
            
            
            if numel(find(intL)) >= numCorrThreshold
                %then add to counter
%                 obs.countL(i,j) = obs.countL(i,j) + 1;
                obs.countL(i,j) = obs.countL(i,j) + numel(find(intL));
            end
            if numel(find(intD)) >= numCorrThreshold
                %then add to counter
%                 obs.countD(i,j) = obs.countD(i,j) +1;
                obs.countD(i,j) = obs.countD(i,j) + numel(find(intD));
            end
            if numel(find(intL)) >= numCorrThreshold && numel(intD) >= numCorrThreshold
                %then add to counter
%                 obs.count(i,j) = obs.count(i,j) + 1;
                obs.count(i,j) = obs.count(i,j) + min([numel(find(intL)) numel(find(intD))]);
            end
        end
    end

    results.obs.corr = corr;
    results.obs.iPairs = iPairs;
    results.obs.jPairs = jPairs;
    results.obs.countn = obs.countn;
    results.obs.countL = obs.countL;
    results.obs.countD = obs.countD;
    results.obs.count = obs.count;
    
    figure; imagesc(results.obs.corr); title('obs.corr'); colorbar
    printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
    figure; imagesc(results.obs.count); title('obs.count'); colorbar
    printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
    figure; imagesc(results.obs.countL); title('obs.countL'); colorbar
    printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
    figure; imagesc(results.obs.countD); title('obs.countD'); colorbar
    printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
    
%---END testing-------------------------------------------------------------

%{
%---setup up data matrix for kmeans clustering----------------
corrValueList = tril(corr,-1);
corrValueList = corrValueList + triu(-1*ones(size(corr))); %get each cell pair's corr value, without the diagonal or repeats by replacing with -1, and logical indexing on next line.
[i j] = find(corrValueList > -1);
idx = sub2ind(size(corrValueList),i,j);
dataframe = corrValueList(idx);
Lcorr = L*L';  %should normalize the spatial values to within sig [0 1] for s, before doing the corr. Else kmeans gets weighted by these higher values too much?
Dcorr = D*D';  %should normalize the spatial values to within sig [0 1] for s, before doing the corr. Else kmeans gets weighted by these higher values too much?
dataframe = [dataframe Lcorr(idx) Dcorr(idx)];
[IDX,C] = kmeans(dataframe,2);
colors = jet(numel(unique(IDX)));
figure; scatter3(dataframe(:,1),dataframe(:,3),dataframe(:,2),50,colors(IDX,:));
%-------------------------------------------------------------
%}

countn = zeros(size(corr));  %setup the results histogram array
countL = zeros(size(corr));
countD = zeros(size(corr));
count = zeros(size(corr));  %setup the results histogram array

h = figure;
ax(1) = subplot(3,1,1);
imagesc(countn); colorbar; title('res.countn'); axis square; drawnow;
ax(2) = subplot(3,1,2);
imagesc(countL); colorbar; title('res.countL'); axis square; drawnow;
ax(3) = subplot(3,1,3);
imagesc(count); colorbar; title('res.count'); axis square; drawnow;
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])

res.s = s;

hbar = waitbar(0,'Please wait...');
stp = max([round(numres)/10 1]);
for t = 1:numres
    %     s = zeros(length(onsets),sz);  %only if you will use the +/- 1frame spike surround procedure (see several lines below)
    res.L = zeros(size(res.s));
    res.D = zeros(size(res.s));
    for c = 1:size(res.s,1)
        res.s(c,:) = res.s(c,randperm(size(res.s,2)));  %shuffle the spike times for each cell
        %         for f = onsets{c}
        %             s(c,:) = s(c,:) + abs((1:sz)-f)<=sig;  %this surrounds the spike with spikes +/-1 frames, for 3 total consecutive frames for ea spike signal.
        %         end
    end
    
    %CHANGE---the following couple lines to a randperm algorithm to make use of the distribution of location and direction values present in the data.
%     res.L(res.s>0) = fix(rand(1,numel(find(res.s>0)))*nLoca)+1;  %randomize the location values for each spike
    for c = 1:size(res.s,1)
        res.L(c,res.s(c,:)>0) = unique(L(c,L(c,:)>0));  %changed to this because randomizing the spike times already randomizes the location of co-occured spikes
    end
    res.D(res.s>0) = fix(rand(1,numel(find(res.s>0)))*nDir)+1; %randomize the direction values for each spike, 4 cardinal directions.
    
    %Smooth spike times with a gaussian (better pvalue calculations)-------
%     res.sthick = gauss_events(res.s,sig);
    
%     %***Smooth location L (+/- 1 frame) and direction D (+/- 1 frame) matrices-----
% %     res.L = L;
% %     res.D = D;
    %***this surrounds the spike with spikes +/-1 frames, for 3 total consecutive frames for ea spike signal.
    for c = 1:size(res.s,1)
        for f = find(res.L(c,:))
            in = find(abs((1:size(res.s,2))-f)<=1);
            res.L(c,in) = res.L(c,f);
            res.D(c,in) = res.D(c,f);
        end
    end
    res.sthick = res.s;
   %} 
% %     figure; imagesc(res.L); title('gauss res.L')
% %     figure; imagesc(res.D); title('gauss res.D')

if t == 1
    figure; 
    ax2(1)=subplot(3,1,1);
    imagesc(res.s); title('res.s')
    ax2(2)=subplot(3,1,2);
    imagesc(res.L); title('res.L')
    ax2(3)=subplot(3,1,3);
    imagesc(res.D); title('res.D')
    linkaxes(ax2); zoom xon;
%     printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
%     return
end
    
    %get the corr value----------------------------------------------------
    res.corr = res.sthick*res.sthick';
%         res.corr = res.s*res.s'; %only if you will test without gaussian spike smoothing
    
    %find temporally correlated cell pairs---------------------------------
    countn = countn + (res.corr>=corr);   
    res.pairs = res.corr>=corr;
    
    %If temporally corr cell pairs are present then find spatial and direction 
    %correlated spike times and set up results arrays----------------------
    [iPairs jPairs] = find(tril(res.pairs,-1));
    if ~isempty(iPairs)
        for c=1:size(iPairs,1)
            %get the corr cell indices i,j
            %         size([i j],1);
            i = iPairs(c);
            j = jPairs(c);
            %find the correlated spike indices for cell i, cell j
            cellIND = res.sthick(i,:) .* res.sthick(j,:);  %element wise multiplication to intersect the spike times for the two corr cells
            %             level = graythresh(cellIND);  %separate the histrogram of corr values into two populations, especially useful if gauss events was used for spike smoothing and there is loads of near-zero values.
            [muhat, sigmahat] = normfit(cellIND);
            level = muhat+(2*sigmahat);
            spkIND = find(cellIND>level);
            
            %CHANGE-- the dfSpk catch algorithm for spikes 1 frame apart in middle---
            dfSpk = find(diff(spkIND) > 1);
            if ~isempty(dfSpk)
%                 spkIND = [spkIND(dfSpk) spkIND(dfSpk+1)];
                spkIND = setdiff(spkIND,spkIND+1);
            end
            %--------------------------------------------------------------
            %test spatial match and direction match at the correlated spike indices in L_res and D_res
            ind.L.i = res.L(i,spkIND);
            ind.L.j = res.L(j,spkIND);
            ind.D.i = res.D(i,spkIND);
            ind.D.j = res.D(j,spkIND);
            %numel n. overlapping spikes
            intL = ismember(ind.L.i(ind.L.i>0),ind.L.j(ind.L.j>0)); %20110901 changed to ismember from intersect, as intersect will not return repeated values
            intD = ismember(ind.D.i(ind.D.i>0),ind.D.j(ind.D.j>0)); %20110901 changed to ismember from intersect, as intersect will not return repeated values
            
%             if numel(find(intL)) >= numCorrThreshold
            if numel(find(intL)) >= obs.countL(i,j)
                %then add to counter
                countL(i,j) = countL(i,j) + 1;
            end
%             if numel(find(intD)) >= numCorrThreshold
            if numel(find(intD)) >= obs.countD(i,j)
                %then add to counter
                countD(i,j) = countD(i,j) + 1;
            end
%             if numel(find(intL)) >= numCorrThreshold && numel(find(intD)) >= numCorrThreshold
            if numel(find(intL)) >= obs.countL(i,j) && numel(find(intD)) >= obs.countD(i,j)
                %then add to counter
                count(i,j) = count(i,j) + 1;
            end
        end
    end
    
%     if mod(t,stp) == 0
%         fprintf('.');
%         set(h,'CurrentAxes',ax(1))
%         imagesc(countn); colorbar; title('res.countn'); axis square; drawnow
%         set(h,'CurrentAxes',ax(2))
%         imagesc(countL); colorbar; title('res.countL'); axis square; drawnow
%         set(h,'CurrentAxes',ax(3))
%         imagesc(count); colorbar; title('res.count'); axis square; drawnow
%     end
    if rem(t,10) == 0
        waitbar(t/numres,hbar); 
    end 
end
close(hbar)
fprintf('\n');
countn
countL
countD
count

% % disp(['countL./countn = ']); countL./countn
% % disp(['countD./countn = ']); countD./countn
p_val = count/numres;
% results = [results p_val(2,1)];
% end
tmp = tril(p_val,-1);
tmp = tmp + triu(ones(size(tmp)));
tmp(nanvalues) = 1;

[i j] = find(tmp < 0.05);
[i j]
size([i j],1)

results.res.countn = countn;
results.res.countL = countL;
results.res.countD = countD;
results.res.count = count;
results.numres = numres;
results.p_val = p_val;
results.pairs = [i j];
results.bigmesh = bigmesh;
results.date=datestr(now,'yyyymmdd-HHMMSS');

figure; imagesc(results.res.countn); title('res.countn'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(results.res.countL); title('res.countL'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(results.res.countD); title('res.countD'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(results.res.count); title('res.count'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(results.p_val); title('results.p-val'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(tmp); title('tmp, res.p-val'); colorbar
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])

[m n] = find(results.obs.count); 
obsPairList = [m n]
% [m n] = ind2sub(size(results.obs.count),idx);
idx = find(results.obs.count);
results.p_val(idx)
signifobsPairList = obsPairList(results.p_val(idx) < 0.05,:)


%contour plot of the ROI pairs---------------------------------------------
figure;
imagesize = size(region.image);
imshow(zeros(imagesize)); colormap(flipud(gray));
hold on;
cl = hsv(length(region.name));
cnt = zeros(1,length(region.contours));
% cnt = zeros(1,length(bigmesh.contours));
roiBinPairList = [];
for c = 1:length(region.contours)
% for c = 1:length(bigmesh.contours)
    cnt(c) = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[1 1 1]);
%     cnt(c) = patch(bigmesh.contours{c}([1:end 1],1),bigmesh.contours{c}([1:end 1],2),[1 1 1]);
    set(cnt(c),'edgecolor',[0.5 0.5 0.5]);
%     set(cnt(c),'edgecolor',[0.5 0.5 0.5]);
%     set(cnt(c),'ButtonDownFcn',['set(txcellnum,''string'',' num2str(c) '); hevButtonDownContours;']);
    for j = 1:length(bigmesh.contours)
        ps = round(region.contours{c});
        binnedROIcoords = round(bigmesh.contours{j});
        [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
        inp = inpolygon(subx,suby,binnedROIcoords(:,1),binnedROIcoords(:,2));
%         if sum(inp(:)) > 0 && ~isempty(intersect(j,results.pairs(:)))  %this one was wrong
        if sum(inp(:)) > 0 && ~isempty(intersect(j,signifobsPairList(:)))  %this one was wrong
            set(cnt(c),'facecolor',[0 0 0]);
            %add c to a counter for results
            roiBinPairList = [roiBinPairList; c j];
        end
    end
end
axis equal
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
box on
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
printfig('png',[fname2base 'contours' datestr(now,'yyyymmdd-HHMMSS') '.png'])

results.roiBinPairList = roiBinPairList;
if ~isfield(region,'userdata') || ~isfield(region.userdata,'corr_pairs_spatialdir')
    region.userdata.corr_pairs_spatialdir{1}=results;
else
    len=length(region.userdata.corr_pairs_spatialdir);
    region.userdata.corr_pairs_spatialdir{len+1}=results;
end
save(fnm,'region')
end

%==========================================================================
function L = setupLocationMatrix(s,bigmesh,region)
locationIndices = unique(bigmesh.location);
centrList={};
centrListAll = [];
%     locations = zeros(size(bigmesh.location));
for i = 1:numel(locationIndices)
    idx = find(bigmesh.location == locationIndices(i));
    locationCode = 1:numel(idx);
    
    %the following code will sort the contours so we can make a matching location code number for each hemisphere
    centrList{i} = zeros(length(idx),2);
    for c = 1:length(idx)
        centrList{i}(c,:)=centroid(bigmesh.contours{idx(c)});   %get centroids of contours
    end
    numContours = idx;   %the existing contour numbers
    centrList{i} = [centrList{i} numContours'];  %attach existing contour numbers to the centroid dataframe
    %     disp(centrList)
    
    [indFirstROI] = find(centrList{i}(:,3) == idx(1));
    if region.orientation.value(1) < centrList{i}(indFirstROI,2)  %determine whether the locationCodes need to be reverse sorted to matchup (depending on whether the region.locations are in same hemipshere or not)
        centrList{i} = sortrows(centrList{i},[1 -2]);   %reverse sort the contour indices so that the locationCode numbers will match up for corresponding locations in each hemisphere.
        centrList{i} = [centrList{i} locationCode'];  %attach the locationCodes
        %     disp(centrList)
    else
        centrList{i} = [centrList{i} locationCode'];
    end
end

for i = 1:length(centrList)
    centrListAll = [centrListAll; centrList{i}];
end

spikeArrayLocations = 1:length(bigmesh.contours);
centrListAll = [centrListAll spikeArrayLocations'];
disp(centrListAll)
L = zeros(size(s));
% [m n] = find(s);
% L(m,n) = m;

for c = 1:size(L,1)
    idx= find(s(c,:));
    [indROI] = find(centrListAll(:,3) == c);
    L(c,idx) = centrListAll(indROI,4);
end
end

%==========================================================================
function D = setupDirectionMatrix(s,bigmesh,region)
locationIndices = unique(bigmesh.location);
wavedir = {};
D = zeros(size(s));
for i = 1:numel(locationIndices)
    locationIndex = locationIndices(i);
    roiIND = find(bigmesh.location == locationIndex);
    
    %wave direction------------------------------------------------------------
    % wavedir = mean(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(~isnan(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians))) * (180/pi);
    centr = round(centroid(region.coords{locationIndex}));
    rowdifference = abs(region.wavedata{locationIndex}.waveprops.wavedirection.origin_coord(1) - region.orientation.value(1));
    if rowdifference < centr(2)
        theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians .* -1; %this flips the coordinate angles in one hemisphere so that the final angles are directly comparable between both hemipheres in the anterior (180deg), medial (270deg), lateral(90deg), and posterior(0deg) directions.
    else
        theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians;
    end
    % wavedir=(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians) .* (180/pi);
%     wavedir{i} = theta .* (180/pi);
    wavedir{i} = theta;
    
    idx = find(wavedir{i} < -2*pi | wavedir{i} > 2*pi);
    if ~isempty(idx)
        error('Some angles outside the interval [-2pi,2pi]')
    end
    
    for wfr = 1:numel(wavedir{i})
        if wavedir{i}(wfr) < pi/4 && wavedir{i}(wfr) >= -pi/4
            wavedir{i}(wfr) = 1;
        elseif wavedir{i}(wfr) < -pi/4 && wavedir{i}(wfr) >= -3*pi/4
            wavedir{i}(wfr) = 2;
        elseif (wavedir{i}(wfr) <= -2*pi && wavedir{i}(wfr) > -3*pi/4) || (wavedir{i}(wfr) <= 2*pi && wavedir{i}(wfr) >= 3*pi/4)
            wavedir{i}(wfr) = 3;
        elseif wavedir{i}(wfr) < 3*pi/4 && wavedir{i}(wfr) >= pi/4
            wavedir{i}(wfr) = 4;
        end
        
        for c = roiIND
            [n] = find(s(c,region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.waveoffsets(wfr)));
            n = n - 1 + region.wavedata{locationIndex}.waveonsets(wfr); 
            D(c,n) = wavedir{i}(wfr);
        end
    end
end


%?(1)make sure wavedirs are on the interval [2*pi, -2*pi] and if not fix angles to be within this angle.
%?(2)make vector with directionCode replacements for 4 cardinal directions for each wave in wavedir
%?(3)for each wave in locationIndex, find [m,n] contours & spiketimes inbetween waveonset:waveoffset. 
%?(4)Using [m,n] idx of spike times, make D array location equal to the theta wavedir{locationIndex}(wfr}
end