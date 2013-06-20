%contourarrayFix.m
%James B. Ackman Mar 2012
%Set up cn array variable equal to the region.contour data that can be passed to calciumdxCentroid, calciumdxDrawCells, and calciumdxManualDelete. Modified from calciumdxFindCells
locationMarkers = unique(region.location);
newRegionLocations = [];
sz = size(region.image);
imageAreaPolygon = [1 1; 1 sz(1); sz(2) sz(1); sz(2) 1];
centr{locationIndex} = [];
areas{locationIndex} = [];
for locationIndex = locationMarkers
        for i = 1:length(cn{locationIndex})
            centr{locationIndex}(i,:) = calciumdxCentroid(cn{locationIndex}{i});
            areas{locationIndex}(i) = polyarea(cn{locationIndex}{i}(:,1),cn{locationIndex}{i}(:,2));
        end
    
    %region.coords is currently the outline of fluoresence for calciumdextran movies (a region).  Formerly, was set to outline of region.image when not doing small regions in the image which is when this contourarraysetup script was made for importing existing contours. Can change script here for temp coords to equal outline of image for polyarea comparison to avoid out of bounds errors.
    in = inpolygon(centr{locationIndex}(:,1),centr{locationIndex}(:,2),imageAreaPolygon(:,1),imageAreaPolygon(:,2));
    %     for c = 1:length(region.name)
    %     for c = 1
    %         if polyarea(region.coords{c}(:,1),region.coords{c}(:,2)) < polyarea(region.coords{num}(:,1),region.coords{num}(:,2))
    %             inoth = inpolygon(centr{locationIndex}(:,1),centr{locationIndex}(:,2),region.coords{craniotomyRegionIndex}(:,1),region.coords{craniotomyRegionIndex}(:,2));
    %             in(find(inoth==1)) = 0;
    %         end
    %     end
    
    f = find(in);
    cntemp = [];
    centrTmp = [];
    areasTmp = [];
    for c = 1:length(f)
        cntemp{c} = cn{locationIndex}{f(c)};
        centrTmp = [centrTmp; centr{locationIndex}(f(c),:)];
        areasTmp = [areasTmp areas{locationIndex}(f(c))];
    end
    cn{locationIndex} = [];
    for c = 1:length(cntemp)
        cn{locationIndex}{c} = cntemp{c};
        newRegionLocations = [newRegionLocations locationIndex];
    end
    centr{locationIndex} = centrTmp;
    areas{locationIndex} = areasTmp;
    
    f = find(region.coords{locationIndex}(:,1) > sz(2));
    region.coords{locationIndex}(f,1) = sz(2);
    f = find(region.coords{locationIndex}(:,1) < 1);
    region.coords{locationIndex}(f,1) = 1;
    f = find(region.coords{locationIndex}(:,2) > sz(1));
    region.coords{locationIndex}(f,2) = sz(1);
    f = find(region.coords{locationIndex}(:,2) < 1);
    region.coords{locationIndex}(f,2) = 1;
    
    
    %     cn{num} = region.contours;
    %     centr{num} = [];
    %     areas{num} = [];
    %     for c = 1:length(cn{num})
    %         centr{num}(c,:) = calciumdxCentroid(cn{num}{c});
    %         areas{num}(c) = polyarea(cn{num}{c}(:,1),cn{num}{c}(:,2));
    %     end
end
region.location = newRegionLocations;
%{
%region.coords is currently the outline of fluoresence for calciumdextran movies (a region).  Formerly, was set to outline of region.image when not doing small regions in the image which is when this contourarraysetup script was made for importing existing contours. Can change script here for temp coords to equal outline of image for polyarea comparison to avoid out of bounds errors.
in = inpolygon(centr{num}(:,1),centr{num}(:,2),region.coords{num}(:,1),region.coords{num}(:,2));
for c = 1:length(region.name)
	if polyarea(region.coords{c}(:,1),region.coords{c}(:,2)) < polyarea(region.coords{num}(:,1),region.coords{num}(:,2))
		inoth = inpolygon(centr{num}(:,1),centr{num}(:,2),region.coords{c}(:,1),region.coords{c}(:,2));
		in(find(inoth==1)) = 0;
	end
end
	
f = find(in);
centr{num} = centr{num}(f,:);
areas{num} = areas{num}(f);
cntemp = [];
for c = 1:length(f)
	cntemp{c} = cn{num}{f(c)};
end
cn{num} = [];
for c = 1:length(cntemp)
	cn{num}{c} = cntemp{c};
end

%}