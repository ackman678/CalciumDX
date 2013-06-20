f = find(spk(num,:)==1);
g = find(dec(num,:)==1);
spk(num,:) = 0;
dec(num,:) = 0;
spk(num,f([1:selev-1 selev+1:end])) = 1;
dec(num,g([1:selev-2 selev:end])) = 1;

hev2PlotTrace