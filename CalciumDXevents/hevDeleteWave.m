function hevDeleteWave
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

disp(['ia = ' num2str(ia)])
disp(['length xdata = ' num2str(length(xdata{2}))])
fidx1 = find(xdata{2}<xdata{2}(ia));
fidx2 = find(xdata{2}>xdata{2}(ia));
% disp(num2str(fidx1))
% disp(num2str(fidx2))

xdata{2} = [xdata{2}(fidx1) xdata{2}(fidx2)];
ydata{2} = [ydata{2}(fidx1) ydata{2}(fidx2)];
set(chi(2),'XData',xdata{2},'YData',ydata{2})

xdata{1} = [xdata{1}(fidx1) xdata{1}(fidx2)];
ydata{1} = [ydata{1}(fidx1) ydata{1}(fidx2)];
set(chi(1),'XData',xdata{1},'YData',ydata{1})

xdata{4} = [xdata{4}(fidx1) xdata{4}(fidx2)];
ydata{4} = [ydata{4}(fidx1) ydata{4}(fidx2)];
set(chi(4),'XData',xdata{4},'YData',ydata{4})
end
