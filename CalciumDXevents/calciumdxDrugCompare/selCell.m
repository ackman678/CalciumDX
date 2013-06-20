function cellFound=selCell(region,imgax)
subplot(imgax)
[xv,yv]=ginput;
cellFound=[];
for i=1:size(region.contours,2);
    in = inpolygon(xv,yv,region.contours{i}(:,1),region.contours{i}(:,2))';
    if length(find(in==1))>0
        cellFound=[cellFound i];
    end
end
