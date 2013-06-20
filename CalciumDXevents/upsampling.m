trNew=[];
for i=1:size(region.traces,1);
    trNew=[trNew; interp(region.traces(i,:),4)];
end
region.traces=trNew;