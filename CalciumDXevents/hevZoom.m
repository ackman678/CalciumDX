% point1 = get(gca,'CurrentPoint');
% finalRect = rbbox;
% point2 = get(gca,'CurrentPoint');
% 
% % idx = get(rtransients,'value');  JBA 8/11/2010
% % region.transients(1,num) = idx;
% 
% 
% candx = sort([point1(1) point2(1)]);
% if range(candx) > 1
%     xlimits = sort([point1(1) point2(1)]);
% else
%     xlimits = [0 size(nt,2)+1];
% end
% hev2PlotTrace

function hevZoom(figHandle,axHandle)
% Allow a line to have its own 'ButtonDownFcn' callback.
% hLine = plot(rand(1,10));
% axHandle = gca;
% set(figHandle,'ButtonDownFcn','disp(''This executes'')');
% set(figHandle,'ButtonDownFcn','hevButtonDown');
set(figHandle,'Tag','DoNotIgnore');
h = zoom;
set(h,'ButtonDownFilter',@mycallback);
set(h,'Enable','on');
setAxesZoomMotion(h,axHandle,'horizontal')
% mouse click on the line
%
function [flag] = mycallback(obj,event_obj)
% If the tag of the object is 'DoNotIgnore', then return true.
objTag = get(obj,'Tag');
if strcmpi(objTag,'DoNotIgnore')
   flag = true;
else
   flag = false;
end
