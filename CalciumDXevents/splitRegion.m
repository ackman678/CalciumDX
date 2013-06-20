function regionsAll = splitRegion(region)
regionsAll = cell(length(region.name),1);
locationMarkers = unique(region.location);
for locationIndex = locationMarkers
   regionsAll{locationIndex}.region = region;
   f = find(region.location == locationIndex);
   regionsAll{locationIndex}.region.traces = zeros(numel(f),size(region.traces,2));
   i = 1:length(f);
   
   regionsAll{locationIndex}.region.contours = {};
   regionsAll{locationIndex}.region.onsets = {};
   regionsAll{locationIndex}.region.offsets = {};
   regionsAll{locationIndex}.region.traces(i,:) = region.traces(f,:);
   if isfield(region,'tracesFilt')
       regionsAll{locationIndex}.region.tracesFilt = zeros(numel(f),size(region.traces,2));
       regionsAll{locationIndex}.region.tracesFilt(i,:) = region.tracesFilt(f,:);
   end
   for i = 1:length(f)
   regionsAll{locationIndex}.region.contours{i} = region.contours{f(i)};
   regionsAll{locationIndex}.region.onsets{i} = region.onsets{f(i)};
   regionsAll{locationIndex}.region.offsets{i} = region.offsets{f(i)};
   end
end