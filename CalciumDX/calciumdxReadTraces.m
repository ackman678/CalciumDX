st = get(dpreaders,'value');
currdir = pwd;
cd(fullfile(calciumdxpath,'TraceReaders'));
[tr, trhalo, param] = feval(mt(st).name(1:end-2),fnm,region,series1);
cd(currdir)

region.traces = tr;
region.halotraces = trhalo;
region.tracereadername = mt(st).name(1:end-2);
region.tracereaderparam = param;

delete(get(fig,'children'));

uicontrol(fig,'Style','text','Units','normalized','String','Signals','Position',[.87 .955 .11 0.03],'FontSize',12,'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);


halo_check = uicontrol(fig,'Style','checkbox','Units','normalized','String','Show halo traces','Position',[.87 .915 .11 0.025],'FontSize',9,...
    'BackgroundColor',[.8 .8 .8],'Callback','calciumdxPlotTrace');
if region.halomode == 0
    set(halo_check,'enable','off');
end

df_check = uicontrol(fig,'Style','checkbox','Units','normalized','String','Calculate DF/F','Position',[.87 .885 .11 0.025],'FontSize',9,...
    'BackgroundColor',[.8 .8 .8],'Callback','calciumdxPlotTrace');

numslider = uicontrol(fig,'Style','slider','Units','normalized','Position',[0.1 0.05 0.74 0.03],'Callback','calciumdxPlotTrace',...
    'Min',1,'Max',length(region.contours),'Sliderstep',[1/length(region.contours) 10/length(region.contours)],'Value',1);

figure(fig)
subplot('position',[0.1 0.6 0.74 0.35]);
hold on;
cl = hsv(length(region.name));
cnt = zeros(1,length(region.contours));
for c = 1:length(region.contours)
    cnt(c) = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
    set(cnt(c),'edgecolor',cl(region.location(c),:));
    set(cnt(c),'ButtonDownFcn',['set(numslider,''value'',' num2str(c) '); calciumdxPlotTrace;']);
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
box on
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);

currdir = pwd;
cd(fullfile(calciumdxpath,'SignalDetectors'));
mt = dir('*.m');
cd(currdir);

st = cell(1,length(mt));
for c = 1:length(mt)
    st{c} = mt(c).name(1:end-2);
    if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'calciumdxSD_')
        st{c} = st{c}(9:end);
    end
end

dummy(1) = uicontrol(fig,'Style','text','Units','normalized','String','Signal detector','Position',[.87 0.8425 .11 0.02],'FontSize',9,...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
dpdetectors = uicontrol(fig,'Style','popupmenu','Units','normalized','String',st,'Position',[.87 .8175 .11 0.025],'FontSize',9,...
    'BackgroundColor',[1 1 1]);
btdetect = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Detect!','Position',[.93 .7725 .05 0.03],'FontSize',9,...
    'Callback','calciumdxDetectSignals');

bnext = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Finish','Position',[.93 .02 .05 .03],'FontSize',9, ...
    'Enable','off','Callback','calciumdxFinish');

region.onsets = cell(1,length(region.contours));
region.offsets = cell(1,length(region.contours));
th = [];
calciumdxPlotTrace;