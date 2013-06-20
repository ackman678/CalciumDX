function region = getActiveFraction(region)
%getActiveFraction.m
%James Ackman, 1/7/2011
%updated to work by location, 5/9/2011
locationMarkers = unique(region.location);
regionsAll = splitRegion(region);
for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    results = getActiveFractionByLocation(tmpregion,locationIndex);
    region.wavedata{locationIndex}.waveprops.waveactvfraction=results;
end

function results = getActiveFractionByLocation(region,locationIndex)
%--------Fraction of ROIs per wave-----------------------------------------
spk = zeros(size(region.traces));
dec = zeros(size(region.traces));
fprintf('\n');
disp(['total no. of ROIs: ' num2str(length(region.contours))])
disp(['no. of active ROIs per wave: '])
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end
actvFraction=[];
for i=1:numel(region.wavedata{locationIndex}.waveonsets)
    %changed the following to count for no. of unique row indices, since each ROI may be active more than once during a wave with the new CCD based recordings
    [row, col] = find(spk(:,region.wavedata{locationIndex}.waveonsets(i):region.wavedata{locationIndex}.waveoffsets(i)) > 0);
%     disp(num2str(unique(row)))
    num = numel(unique(row));
    disp(num2str(num));
    tmp = num / length(region.contours);
    actvFraction=[actvFraction tmp];
end
disp(['Active fraction: ' num2str(actvFraction)])
results=actvFraction;