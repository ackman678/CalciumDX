%Repeat delete of ROIs
butt = 1;  %left mouse click is '1'
while butt < 2 %right mouse click is '2'
    [x,y,butt] = ginput(1);
    if butt > 1
        return
    end
    
    %Get matrix with the x, y centroid coords for each ROI together with its 'region.location' index and the ROI no. for that region
    matr = [];
    for c = 1:length(cn)
        for d = 1:length(cn{c})
            if polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) > lowar(c) & polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) < highar(c)
                matr = [matr; [calciumdxCentroid(cn{c}{d}) c d]];
            end
        end
    end
    
    %Make a min distance calculation for all ROIs against the x,y coords of the ROI to be deleted
    dst = sum((matr(:,1:2)-repmat([x y],size(matr,1),1)).^2,2);
    [dummy i] = min(dst);
    
    %Delete the contour array for the ROI to be deleted and setup temporary cn variable
    tempcn = [];
    for d = [1:matr(i,4)-1 matr(i,4)+1:length(cn{matr(i,3)})]
        tempcn{length(tempcn)+1} = cn{matr(i,3)}{d};
    end
    
    %Delete the value from region.location if it exists yet (for isimported=TRUE)
    if isfield(region,'location')
        newRegionLocations = [];
        for d = [1:i-1 i+1:size(matr,1)]
            newRegionLocations = [newRegionLocations matr(d,3)];
        end
        region.location = newRegionLocations;
    end
    
    %Finally save the set of new ROI contours to the cn variable
    cn{matr(i,3)} = [];
    for d = 1:length(tempcn)
        cn{matr(i,3)}{d} = tempcn{d};
    end
    
    tmp = num;
    num = matr(i,3);
    calciumdxDrawCells
    num = tmp;
end