function hevMoveWaveOnset
[x y butt] = ginput(1);
if butt > 1
    return
end
x = round(x);
y = round(y);
chi=get(gca,'Children');
xdata=get(chi,'XData');
ydata=get(chi,'YData');

set1 = sort([xdata{1} xdata{4}]);
fidx1 = find(set1<x);
fidx2 = find(set1>x);
% disp(fidx1)
% disp(val)

if isempty(fidx1)
    fidx1 = 1;
else
    fidx1 = max([1 set1(fidx1(end))]);
end
if isempty(fidx2)
    fidx2 = length(xdata{5});
else
    fidx2 = min([length(xdata{5}) set1(fidx2(1))]);
end

finterval = fidx1:fidx2;
[c,ia,ib] = intersect(xdata{2},finterval);
xdata{2}(ia) = x;
ydata{2}(ia) = ydata{5}(x);
set(chi(2),'XData',xdata{2},'YData',ydata{2})

% %algo for adding new onsets
% ydata=get(chi,'YData');
% ydata{2} = ([ydata{2} ydata{5}(x)]);
% xdata{2} = sort([xdata{2} x]);
% set(chi(2),'XData',xdata{2},'YData',ydata{2});


% disp(finterval)

% if x < 1
%     return
% end
% 
% f = find(spk(num,:)==1);
% g = find(dec(num,:)==1);
% 
% if selev > 1 & x <= g(selev-1)
%     return
% end
% if x >= g(selev)
%     return
% end
% 
% spk(num,f(selev)) = 0;
% spk(num,x) = 1;
% hev2PlotTrace
end
%}