for c = 1:length(reg)
    region.name{c} = get(inpt(c),'String');
end

delete(txlab)
delete(inpt)
delete(bnext)

handl = {};

islocal = 0;

det_loc = uicontrol('Style','pushbutton','Units','normalized','String','Filter','Position',[.87 .535 .05 .03],'FontSize',9, ...
    'Callback','calciumdxFilter');

det_view = uicontrol('Style','pushbutton','Units','normalized','String','View','Position',[.93 .535 .05 .03],'FontSize',9, ...
    'Callback','calciumdxViewLoc','enable','off');

num = 1;
loca = region.image;
region.filtername = 'none';
region.filterparam = [];

calciumdxInputParams;
