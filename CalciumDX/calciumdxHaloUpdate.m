delete(halo_hands);
halo_hands = [];
halos = cell(1,length(region.contours));
region.haloarea = str2num(get(inpthaloar,'string'));
for c = 1:length(region.contours);
    ct = repmat(centroid(region.contours{c}),size(region.contours{c},1),1);
    halos{c} = (region.contours{c}-ct)*sqrt((1+region.haloarea))+ct;
    halo_hands(c) = plot(halos{c}([1:end 1],1),halos{c}([1:end 1],2),'color',cl(region.location(c),:),'LineWidth',1,'LineStyle',':');
end
zoom on