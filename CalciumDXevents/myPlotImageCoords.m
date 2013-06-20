%myPlotImageCoords
%JBA, Friday, January 20, 2012 4:51 PM

%{
%plot multiple images with overlays for each---------
figure;
ax(1) = subplot(3,1,1)
imshow(mat2gray(region.image))
% colormap(ax(1),gray)
hold on
for numcoords = 1:length(region.coords)
    plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
end

ax(2) = subplot(3,1,2)
imshow(mat2gray(region.image))
hold on
for c = 1:size(region.contours,2)
   plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0 1 1]);
end

ax(3) = subplot(3,1,3)
[szX,szY] = size(region.image);
Adummy=zeros(szX,szY);
A2=mat2gray(Adummy);
imshow(A2); 
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
for numcoords = 1:length(region.coords)
    plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
end

%}



%just plot one image with overlay for both-----------
figure;
imshow(mat2gray(region.image))
% colormap(ax(1),gray)
hold on
for numcoords = 1:length(region.coords)
    plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
end

for c = 1:size(region.contours,2)
   plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0 1 1]);
end