%calciumdx help 20100904
%avg image
% imgax = subplot('position',[0.02 0.02 0.82 0.96]);

hold off
imagesc(a);
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap hot
% colormap gray

%loca
%calciumdx help 20100904
% imgax = subplot('position',[0.02 0.02 0.82 0.96]);

hold off
imagesc(loca);
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap hot


%Sort ROIs in ascending order by approximately the time in which the pixels were scanned by the laser.
% function centroid = calciumdxCentroid(coords)
centroid=zeros(numel(region.contours),2);
for i=1:numel(region.contours)
coords=region.contours{i};

if prod(size(coords))==0
   cx = NaN;
   cy = NaN;
else
   m = [coords; coords(1,:)];
   x = m(:,1);
   y = m(:,2);
   
   a = (sum(x(1:end-1).*y(2:end)) - sum(x(2:end).*y(1:end-1)))/2;
   cx = sum((x(1:end-1)+x(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*a);
   cy = sum((y(1:end-1)+y(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*a);
end

centroid(i,:) = [cx cy];
end

centroid(:,2)=centroid(:,2).^2;
[B idx]=sortrows(centroid,2);

for j=1:numel(idx)
    i=idx(j);
    disp(i)
    %do something
end
