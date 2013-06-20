for c = 1:length(reg)
    region.name{c} = get(inpt(c),'String');
end

delete(txlab)
delete(inpt)
delete(bnext)

handl = {};

islocal = 0;

det_tx1 = uicontrol('Style','text','Units','normalized','String','Image filter','Position',[.87 .60 .11 0.02],'FontSize',9,...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);

currdir = pwd;
cd(fullfile(calciumdxpath,'ImageFilters'));
mt = dir('*.m');
cd(currdir);

st = cell(1,length(mt));
for c = 1:length(mt)
    st{c} = mt(c).name(1:end-2);
    if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'calciumdxIF_')
        st{c} = st{c}(9:end);
    end
end

dpfilters = uicontrol('Style','popupmenu','Units','normalized','String',st,'Position',[.87 .5725 .11 0.025],'FontSize',9,...
    'BackgroundColor',[1 1 1]);

det_loc = uicontrol('Style','pushbutton','Units','normalized','String','Filter','Position',[.87 .535 .05 .03],'FontSize',9, ...
    'Callback','calciumdxFilter');

det_view = uicontrol('Style','pushbutton','Units','normalized','String','View','Position',[.93 .535 .05 .03],'FontSize',9, ...
    'Callback','calciumdxViewLoc','enable','off');