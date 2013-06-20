rtrans = get(rtransients,'value');
if rtrans == 1
    region.transients(1,num) = 1;
elseif rtrans == 2
    region.transients(1,num) = 2;
elseif rtrans == 3
    region.transients(1,num) = 3;
elseif rtrans == 4
    region.transients(1,num) = 4;
elseif rtrans == 5
    region.transients(1,num) = 5;
end
    

% if get(radio(1),'value') == 1
% region.transients(1,num) = 1
% end
% if get(radio(2),'value') == 1
% region.transients(1,num) = 2
% end
% if get(radio(3),'value') == 1
% region.transients(1,num) = 3
% end
% if get(radio(4),'value') == 1
% region.transients(1,num) = 4
% end

hev2PlotTrace;
