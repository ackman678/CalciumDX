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

onsets = cell(1,length(region.contours));
offsets = cell(1,length(region.contours));
param = [];

region.onsets = onsets;
region.offsets = offsets;
region.detectorname = '';
region.detectorparam = param;

bnext = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Finish','Position',[.93 .02 .05 .03],'FontSize',9, ...
    'Callback','calciumdxFinish');

region.onsets = cell(1,length(region.contours));
region.offsets = cell(1,length(region.contours));
th = [];
calciumdxPlotTrace;