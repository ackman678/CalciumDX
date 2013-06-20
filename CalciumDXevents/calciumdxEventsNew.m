clear;

opengl neverselect;
sz = get(0,'screensize');
fig = figure('Name','calciumdxEvents','NumberTitle','off','MenuBar','none','position',[1 0.15*sz(4) sz(3) 0.7*sz(4)],'doublebuffer','on');

[filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open')

if ~isstr(filename)
    delete(gcf)
    clear
end
fnm = [pathname filename];
load(fnm)

tr = region.traces;
nt = [];
for c = 1:size(tr,1)
    nt(c,:) = dfoverf(tr(c,:))*100;
end

spk = zeros(size(nt));
dec = zeros(size(nt));
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end

xlimits = [0 size(nt,2)+1];

trax = subplot('position',[0.03 0.20 0.94 0.75]);
box on
set(gca,'buttondownfcn','hevZoom')

uicontrol('Style','text','Units','normalized','String','Cell #','Position',[.05 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
txcellnum = uicontrol('Style','edit','Units','normalized','String','1','Position',[.11 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','hevPlotTrace');
bgoto = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.17 .05 .05 0.04],'FontSize',12,...
    'Callback','hevPlotTrace');

progtx = uicontrol('Style','text','Units','normalized','String','','Position',[.70 .05 .25 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);

currdir = pwd;
cd('C:\Programmi\MATLAB\R2006b\work\calciumdx\SignalDetectors');
mt = dir('*.m');
cd(currdir);

st = cell(1,length(mt));
for c = 1:length(mt)
    st{c} = mt(c).name(1:end-2);
    if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'calciumdxEvent_')
        st{c} = st{c}(9:end);
    end
end

dummy(2) = uicontrol('Style','text','Units','normalized','String','Trace reader','Position',[.80 0.12 .11 0.04],'FontSize',9,...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
dpreaders = uicontrol('Style','popupmenu','Units','normalized','String',st,'Position',[.80 .1 .11 0.025],'FontSize',9,...
    'BackgroundColor',[1 1 1]);

bdetect1 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current usual','Position',[.25 .05 .14 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxdettrial(tr(num,:)); spk(num,s) = 1; dec(num,d) = 1; hevPlotTrace;');
bdetect2 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Filt','Position',[.40 .1 .14 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP(region,num,''no''); spk(num,s) = 1; dec(num,d) = 1; hevPlotTrace;');
bdetect3 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Filt Noise','Position',[.40 .005 .14 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP(region,num,''yes''); spk(num,s) = 1; dec(num,d) = 1; hevPlotTrace;');
bdetect4 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Der','Position',[.61 .1 .14 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHPDer(region,num,''no''); spk(num,s) = 1; dec(num,d) = 1; hevPlotTrace;');
bdetect5 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Der Noise','Position',[.61 .005 .14 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHPDer(region,num,''yes''); spk(num,s) = 1; dec(num,d) = 1; hevPlotTrace;');
bdetect = uicontrol('Style','pushbutton','Units','normalized','String','Detect all','Position',[.80 .02 .07 0.04],'FontSize',12,...
    'Callback','calciumdxDispTrialNew');
bdeleteall = uicontrol('Style','pushbutton','Units','normalized','String','Delete events','Position',[.54 .05 .07 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; hevPlotTrace;');
bsave = uicontrol('Style','pushbutton','Units','normalized','String','Save','Position',[.90 .02 .07 0.04],'FontSize',12,...
    'Callback','hevSave');

hevPlotTrace

