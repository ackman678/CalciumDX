function centroid = calciumdxCentroid(coords)

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

centroid = [cx cy];