%hevFuzzyDeleteAll
f = find(spk(num,:)==1);
fr=f(selev);
sz=size(region.traces);
for i=1:length(region.contours)
    fr_range=max([1 fr-3]):min([sz(2) fr+3]);
    f = find(spk(i,:)==1);
    g = find(dec(i,:)==1);
    for v=fr_range
        if any(f == v)
            spk(i,v) = 0;
            [dummy v2]=find(f == v);
            dec(i,g(v2)) = 0;
%             fr=v;
        end
    end
    
    if sum(spk(i,:)) == 0
        region.transients(1,i) = 1;
    end
end
hevPlotTrace;