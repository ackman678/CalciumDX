[x y butt] = ginput(1);
if butt > 1
    return
end
x = round(x);

if x > size(nt,2)
    return
end

f = find(spk(num,:)==1);
g = find(dec(num,:)==1);

if selev < length(f) & x >= f(selev+1)
    return
end
if x <= f(selev)
    return
end

dec(num,g(selev)) = 0;
dec(num,x) = 1;
hevPlotTrace