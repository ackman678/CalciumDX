set(gcf,'CurrentAxes',imgax1); dx1 = ginput(3);
set(gcf,'CurrentAxes',imgax2); dx2 = ginput(3);
dxy = mean(dx2 - dx1);
dx = dxy(1);
dy = dxy(2);

region.contours = cell(1,length(region1.contours));
for c = 1:length(region1.contours); 
    region.contours{c}(:,1) = region1.contours{c}(:,1) - dx; 
    region.contours{c}(:,2) = region1.contours{c}(:,2) - dy; 
end; 

allH=[];
for i=1:length(region.contours);
    hold on;
    h=plot(region.contours{i}([1:end 1],1),region.contours{i}([1:end 1],2),'r-');
    allH=[allH h];
end