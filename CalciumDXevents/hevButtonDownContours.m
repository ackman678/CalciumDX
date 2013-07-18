% num = str;
% % num = 100;
% set(txcellnum,'string',num2str(num));

num = str2double(get(txcellnum,'string'));
set(numslider,'Value',num);
figure(fig);
set(fig,'CurrentAxes',ax(3));
hevPlotTrace;
