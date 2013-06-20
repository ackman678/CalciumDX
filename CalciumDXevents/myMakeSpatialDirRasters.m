function [s, L, D, sthick, Lthick, Dthick, ntSpk,sOnsets] = myMakeSpatialDirRasters(region,fnm,locationIndexPair)
%[s, L, D, sthick, Lthick, Dthick, ntSpk,sOnsets] = myMakeSpatialDirRasters(region,fname)
%2011-12-21 James B. Ackman
if nargin < 3 || isempty(locationIndexPair), locationIndexPair = [2 3]; end

% load matlab3.mat
% load corrspatialtest_111108.mat

fname2base = [fnm(1:end-4) 'dCorrTDL_'];

sig = 1;
nDir = 4;

%---Create spike array movie---
A = myMakeContourMovieWavesOnsets([],region);

%---Dilate ROIs to connect neighbors and exclude isolated pixels having no neighbors nearby (activations not within the activity wave front)
strel_sz = getStrelsize(region);
[A2,A3] = dilateContours(A,strel_sz);

%---Setup new spatially binned ROI spike array containing spikes from array found from the dilated spike array A2
[bigmesh] = makeMatchingROImeshgrid(region,locationIndexPair);
s = intersectSpikeArrayContours(bigmesh.contours,A2);

sthick=zeros(size(s));
for c = 1:size(sthick,1)
    for findSpks = find(sthick(c,:))
        in = find(abs((1:size(sthick,2))-findSpks)<=1);  %smooth the spikes with +/- 1 signal
        sthick(c,in) = sthick(c,findSpks);
    end
end

ntSpk=zeros(size(s));
locationIndex = 2;
for c = 1:100
    for findSpks = find(s(c,:))
        ind = find(region.wavedata{locationIndex}.waveonsets <= findSpks);
        if ~isempty(ind)
            if findSpks >= region.wavedata{locationIndex}.waveonsets(ind(end)) && findSpks <= region.wavedata{locationIndex}.waveoffsets(ind(end))
                ntSpk(c,findSpks) = abs(findSpks-region.wavedata{locationIndex}.waveoffsets(ind(end)))/abs(region.wavedata{locationIndex}.waveonsets(ind(end))-region.wavedata{locationIndex}.waveoffsets(ind(end)));
                in = find(abs((1:size(s,2))-findSpks)<=1);  %smooth the spikes with +/- 1 signal
                s(c,in) = s(c,findSpks);
            end
        end
    end
end

locationIndex = 3;
for c = 101:200
    for findSpks = find(s(c,:))
        ind = find(region.wavedata{locationIndex}.waveonsets <= findSpks);
        if ~isempty(ind)
            if findSpks >= region.wavedata{locationIndex}.waveonsets(ind(end)) && findSpks <= region.wavedata{locationIndex}.waveoffsets(ind(end))
                ntSpk(c,findSpks) = abs(findSpks-region.wavedata{locationIndex}.waveoffsets(ind(end)))/abs(region.wavedata{locationIndex}.waveonsets(ind(end))-region.wavedata{locationIndex}.waveoffsets(ind(end)));
                in = find(abs((1:size(s,2))-findSpks)<=1);  %smooth the spikes with +/- 1 signal
                s(c,in) = s(c,findSpks);
            end
        end
    end
end

sOnsets = zeros(2,size(s,2));
locationIndices = [2 3];
for i = 1:length(locationIndices)
    waveonsets = region.wavedata{locationIndices(i)}.waveonsets;
    sOnsets(i,waveonsets) = 1;
end

for c = 1:size(sOnsets,1)
    for findSpks = find(sOnsets(c,:))
        in = find(abs((1:size(sOnsets,2))-findSpks)<=1);  %smooth the spikes with +/- 1 frame signal
        sOnsets(c,in) = sOnsets(c,findSpks);
    end
end

% figure; imagesc(sOnsets)
% 
% myRasterPlot(sOnsets)


figure; imagesc(s); title('s'); colorbar
% printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
%--------------------------------------------------------------------------------

%--Setup location spike matrix----------------
%L = region.locations  %need to setup
L = setupLocationMatrix(s,bigmesh,region);
Lthick = setupLocationMatrix(sthick,bigmesh,region);


%--Setup direction spike matrix----------------
%D = region.directions  %need to setup
D = setupDirectionMatrix(s,bigmesh,region);
Dthick = setupDirectionMatrix(sthick,bigmesh,region);

figure; imagesc(L); title('L'); colorbar
% printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
figure; imagesc(D); title('D'); colorbar
% printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])


figure;
ax(1)= subplot(3,1,1);
myRasterPlot(sOnsets)
ax(2) = subplot(3,1,2); imagesc(L); title('L'); colorbar
ax(3) = subplot(3,1,3); imagesc(D); title('D'); colorbar
linkaxes(ax,'x')
zoom xon
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
    
    idx = find(wavedir{i} < -pi | wavedir{i} > pi);
    if ~isempty(idx)
        error('Some angles outside the interval [-pi,pi]')
    end
    
    for wfr = 1:numel(wavedir{i})
%Uncomment following lines to bin directions into 1 of 4 cardinal directions
%         if wavedir{i}(wfr) < pi/4 && wavedir{i}(wfr) >= -pi/4
%             wavedir{i}(wfr) = 1;
%         elseif wavedir{i}(wfr) < -pi/4 && wavedir{i}(wfr) >= -3*pi/4
%             wavedir{i}(wfr) = 2;
%         elseif (wavedir{i}(wfr) <= -2*pi && wavedir{i}(wfr) > -3*pi/4) || (wavedir{i}(wfr) <= 2*pi && wavedir{i}(wfr) >= 3*pi/4)
%             wavedir{i}(wfr) = 3;
%         elseif wavedir{i}(wfr) < 3*pi/4 && wavedir{i}(wfr) >= pi/4
%             wavedir{i}(wfr) = 4;
%         end
        
        for c = roiIND
            [n] = find(s(c,region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.waveoffsets(wfr)));
            n = n - 1 + region.wavedata{locationIndex}.waveonsets(wfr);
            D(c,n) = 2*pi + wavedir{i}(wfr);
        end
    end
end


%?(1)make sure wavedirs are on the interval [pi, -pi] and if not fix angles to be within this angle.
%?(2)make vector with directionCode replacements for 4 cardinal directions for each wave in wavedir
%?(3)for each wave in locationIndex, find [m,n] contours & spiketimes inbetween waveonset:waveoffset.
%?(4)Using [m,n] idx of spike times, make D array location equal to the theta wavedir{locationIndex}(wfr}
end