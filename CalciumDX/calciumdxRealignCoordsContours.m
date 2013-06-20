%calciumdxRealignCoordsContours
%James B. Ackman, 2011-2012
[x,y] = ginput(2);
dx = x(1) - x(2);
dy = y(1) - y(2);

if isimported > 0
    locationMarkers = unique(region.location);
    for num = locationMarkers
        for c = 1:length(cn{num})
            cn{num}{c}(:,1) = cn{num}{c}(:,1) - dx;
            cn{num}{c}(:,2) = cn{num}{c}(:,2) - dy;
        end
        region.coords{num}(:,1) = region.coords{num}(:,1) - dx;
        region.coords{num}(:,2) = region.coords{num}(:,2) - dy;
    end
    
    contourarrayFix;
    
%     for c = 1:length(region.contours)
%         region.contours{c}(:,1) = region.contours{c}(:,1) - dx;
%         region.contours{c}(:,2) = region.contours{c}(:,2) - dy;
%     end
    
else
    for c = 1:length(cn{num})
        cn{num}{c}(:,1) = cn{num}{c}(:,1) - dx;
        cn{num}{c}(:,2) = cn{num}{c}(:,2) - dy;
    end
    for numcoords = 1:length(region.coords)
        if prod(max(region.coords{numcoords})) ~= prod(size(region.image))
            region.coords{numcoords}(:,1) = region.coords{numcoords}(:,1) - dx;
            region.coords{numcoords}(:,2) = region.coords{numcoords}(:,2) - dy;
        end
    end
end
calciumdxDrawCells;