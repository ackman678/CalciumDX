%contourarraysetup2
%James B. Ackman, 2008-2012
%Set up cn array variable equal to the region.contour data that can be passed to calciumdxCentroid, calciumdxDrawCells, and calciumdxManualDelete. Modified from calciumdxFindCells
locationMarkers = unique(region.location);
newRegionLocations = [];
sz = size(region.image);
imageAreaPolygon = [1 1; 1 sz(1); sz(2) sz(1); sz(2) 1];

for locationIndex = locationMarkers
    fInd = find(region.location == locationIndex);
    for i = 1:length(fInd)
        cn{locationIndex}{i} = region.contours{fInd(i)};
        centr{locationIndex}(i,:) = calciumdxCentroid(cn{locationIndex}{i});
        areas{locationIndex}(i) = polyarea(cn{locationIndex}{i}(:,1),cn{locationIndex}{i}(:,2));
    end
%     for c = 1:length(cntemp)
%         cn{locationIndex}{c} = cntemp{c};
%         newRegionLocations = [newRegionLocations locationIndex];
%     end 

in = inpolygon(centr{locationIndex}(:,1),centr{locationIndex}(:,2),imageAreaPolygon(:,1),imageAreaPolygon(:,2));
%     for c = 1:length(region.name)
%     for c = 1
%         if polyarea(region.coords{c}(:,1),region.coords{c}(:,2)) < polyarea(region.coords{num}(:,1),region.coords{num}(:,2))
%             inoth = inpolygon(centr{locationIndex}(:,1),centr{locationIndex}(:,2),region.coords{craniotomyRegionIndex}(:,1),region.coords{craniotomyRegionIndex}(:,2));
%             in(find(inoth==1)) = 0;
%         end
%     end

f = find(in);



f = find(region.coords{locationIndex}(:,1) > sz(2));
coords{locationIndex}(f,1) = sz(2);
f = find(region.coords{locationIndex}(:,1) < 1);
coords{locationIndex}(f,1) = 1;
f = find(region.coords{locationIndex}(:,2) > sz(1));
coords{locationIndex}(f,2) = sz(1);
f = find(region.coords{locationIndex}(:,2) < 1);
coords{locationIndex}(f,2) = 1;

end