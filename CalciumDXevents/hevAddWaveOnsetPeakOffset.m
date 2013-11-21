function hevAddWaveOnsetPeakOffset
[xAll yAll butt] = ginput(3);
if butt > 1
    return
end

xAll = round(xAll);
yAll = round(yAll);
chi=get(gca,'Children');
xdata=get(chi,'XData');
ydata=get(chi,'YData');

len = size(xdata{end},2);

%algo for adding new onsets
x=xAll(1);
if x < 1
	x = 1;
elseif x > len
	x = len;
else
	x = x;
end
fidx1 = find(xdata{2}<x);
fidx2 = find(xdata{2}>x);
xdata{2} = sort([xdata{2} x]);
ydata{2} = ([ydata{2}(fidx1) ydata{5}(x) ydata{2}(fidx2)]);
set(chi(2),'XData',xdata{2},'YData',ydata{2});

%algo for adding new peaks
x=xAll(2);
if x < 1
	x = 1;
elseif x > len
	x = len;
else
	x = x;
end
fidx1 = find(xdata{4}<x);
fidx2 = find(xdata{4}>x);
xdata{4} = sort([xdata{4} x]);
ydata{4} = ([ydata{4}(fidx1) ydata{5}(x) ydata{4}(fidx2)]);
set(chi(4),'XData',xdata{4},'YData',ydata{4});

%algo for adding new offsets
x=xAll(3);
if x < 1
	x = 1;
elseif x > len
	x = len;
else
	x = x;
end
fidx1 = find(xdata{1}<x);
fidx2 = find(xdata{1}>x);
xdata{1} = sort([xdata{1} x]);
ydata{1} = ([ydata{1}(fidx1) ydata{5}(x) ydata{1}(fidx2)]);
set(chi(1),'XData',xdata{1},'YData',ydata{1});
