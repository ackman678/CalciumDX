f = find(spk(num,:)==1);
g = find(dec(num,:)==1);
spk(num,f(selev)) = 0;
dec(num,g(selev)) = 0;

if sum(spk(num,:)) == 0
    region.transients(1,num) = 1;
end
hev2PlotTrace;