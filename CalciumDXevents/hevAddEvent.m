tic
[x,y,butt] = ginput(1);
toc
if butt > 1
    return
end
x = round(x);
if x < 1 || x > size(nt,2)
    return
end
if sum(spk(num,1:x-1)) > sum(dec(num,1:x-1))
    return
end
% if length(find(spk{num}<x)) > length(find(dec{num}<x))
%     return
% end


[x2,y,butt] = ginput(1);
if butt > 1
    return
end
x2 = round(x2);
if x2 < 1 || x2 > size(nt,2)
    return
end
if sum(spk(num,1:x2-1)) > sum(dec(num,1:x2-1))
    return
end
if x2 <= x
    return
end

spk(num,x) = 1;
dec(num,x2) = 1;
% spk{num} = sort([spk{num} x]);
% dec{num} = sort([dec{num} x2]);
if region.transients(1,num) == 1;
    region.transients(1,num) = 4;
end

% delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6); delete(h7); delete(h8); delete(h9); delete(h10);
% clear trmenu cmenu coffmenu h1 h2 h3 h4 h5 h6 h7 h8 h9 h10
hev2PlotTrace
