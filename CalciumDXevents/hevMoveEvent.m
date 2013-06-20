[x y butt] = ginput(1);
if butt > 1
    return
end
x = round(x);

if x < 1
    return
end

f = find(spk(num,:)==1);
g = find(dec(num,:)==1);

if selev > 1 & x <= g(selev-1)
    return
end
if x >= g(selev)
    return
end

spk(num,f(selev)) = 0;
spk(num,x) = 1;
hev2PlotTrace
