%updated by James Ackman 2010-2012.
%Updated by James Ackman 22/08/07.
%hevPlotTrace updated 05/04/07 JBA
%clear;

%TODO: change to function based code
%TODO: update, add flexibility for using external plugins for detection routines (like ICA/PCA)
%TODO: add user preferences for gui
%TODO: add user preferences for detection routine argin (like filter order, pass band, window)

%If detecting a positive signal (like with OGB1AM, Fluo4AM, etc) set the
%following line to '-1' to invert the raw traces just for detection purposes. Otherwise leave at '+1' since the autodetect
%scripts are written to detect a decreasing signal by default.

trSign = -1;

%opengl neverselect;
sz = get(0,'screensize');
% fig = figure('Name','calciumdxEvents','NumberTitle','off','MenuBar','none','position',[1 0.15*sz(4) sz(3) 0.7*sz(4)],'doublebuffer','on');
fig = figure('Name','calciumdxEvents','NumberTitle','off','MenuBar','none','position',[1 sz(4) sz(3) sz(4)],'doublebuffer','on');

if ~exist('region','var')
    try
        load calciumdxprefs
    end
    if exist('pathname','var')
        [filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open',pathname);
        if ~isstr(filename)
            delete(gcf)
            clear
        end
        save('calciumdxprefs.mat', 'pathname')
    else
        [filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open');
        if ~isstr(filename)
            delete(gcf)
            clear
        end
        save('calciumdxprefs.mat', 'pathname')
    end
    
    fnm = [pathname filename];
    load(fnm)
end

tr = region.traces;
nt = zeros(size(tr));
for c = 1:size(tr,1)
    nt(c,:) = dfoverf(tr(c,:))*100;
end

maxY = max(nt(:));
minY = min(nt(:));
ylimits = [minY maxY];

spk = zeros(size(nt));
dec = zeros(size(nt));
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end

% trNew=[];
% for i=1:size(region.traces,1);
%     trNew=[trNew; interp(region.traces(i,:),4)];  %interp trace for high frequency transient detection algorithms.
% end


%The following will determine if region.transients exists, and if not will create it
tmp = isfield(region,'transients');
if tmp == 0
region.transients = ones(1,size(tr,1));
end

xlimits = [0 size(nt,2)+1];



%Subplot of contour data-----------------------------------------------------
% ax(1)=subplot('position',[0.03 0.81 0.95 0.18]);
fig2 = figure;
imagesc(region.image); colormap(gray)
hold on;
cl = hsv(length(region.name));
cnt = zeros(1,length(region.contours));
for c = 1:length(region.contours)
    cnt(c) = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
    set(cnt(c),'edgecolor',cl(region.location(c),:));
%     set(cnt(c),'ButtonDownFcn','str = num2str(c); hevButtonDownContours;');
    set(cnt(c),'ButtonDownFcn',['set(txcellnum,''string'',' num2str(c) '); hevButtonDownContours;']);
%     set(cnt(c),'ButtonDownFcn',['set(numslider,''value'',' num2str(c) '); calciumdxPlotTrace;']);
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
box on
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);


% colormap plot of cell fluorescence over time---------------------
figure(fig);
ax(2)=subplot('position',[0.03 0.70 0.94 0.10]);
imagesc(nt); colormap(jet)
ylabel('Cell #')
% set(gca,'ytick',[1 (1:fix(length(region.contours)/100))*100])
set(gca,'ydir','reverse')
% set(gca,'xtick',[1 size(nt,2)]);
% colorbar('North')
% mins = fix(size(nt,2)*tmpres/60);
% secs = fix(size(nt,2)*tmpres - fix(size(nt,2)*tmpres/60)*60);
% if secs == 0
%     str = '00';
% elseif secs < 10
%     str = ['0' num2str(secs)];
% else
%     str = num2str(secs);
% end
% %box on
% set(gca,'xticklabel',{'0', [num2str(mins) ':' str]});
% xlabel('Time (min)')






%Subplot of mean trace for all cells-----------------------------------
ax(3)=subplot('position',[0.03 0.61 0.94 0.08]);
plot(mean(nt));
xlim(xlimits);
set(gca,'xtick',[]);
%set(gca,'ytick',[]);

trax = subplot('position',[0.03 0.15 0.94 0.45]);
box on
% set(gca,'buttondownfcn','hevZoom')
% set(fig,'KeyPressFcn','hevButtonDown')

numslider = uicontrol('Style','slider','Units','normalized','Position',[0.1 0.1 0.74 0.03],'Callback','hevPlotTrace',...
    'Min',1,'Max',length(region.contours),'Sliderstep',[1/length(region.contours) 10/length(region.contours)],'Value',1);



uicontrol('Style','text','Units','normalized','String','Cell #','Position',[.05 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
txcellnum = uicontrol('Style','edit','Units','normalized','String','1','Position',[.11 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','hevPlotTrace');
bgoto = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.17 .05 .05 0.04],'FontSize',12,...
    'Callback','hevPlotTrace');

% progtx = uicontrol('Style','text','Units','normalized','String','','Position',[.70 .05 .25 0.04],'FontSize',12,'FontWeight','Bold',...
%     'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);

%Row1 buttons
%popupmenu replacement for detection routine selection, JBA 07/15/2010
popupDetectList= cell(1,5);
popupDetectList{1} = 'calciumdxdettrial'; 
popupDetectList{2} = 'calciumdxdettrialWaves';
popupDetectList{3} = 'calciumdxEvent_DetSingTrHannFast';
popupDetectList{4} = 'calciumdxEvent_DetSingTrHP';
popupDetectList{5} = 'calciumdxEvent_DetSingTrHP_HF';
popupDetect = uicontrol('Style','popupmenu','Units','normalized','String',popupDetectList,'Position',[.28 .05 .11 0.03],'FontSize',9,...
    'BackgroundColor',[1 1 1]);

bdetect = uicontrol('Style','pushbutton','Units','normalized','String','Detect current','Position',[.40 .05 .07 0.03],'FontSize',9,...
    'Callback','calciumdxeventDetect');
    
bdetect1 = uicontrol('Style','pushbutton','Units','normalized','String','Detect all','Position',[.50 .05 .07 0.03],'FontSize',9,...
    'Callback','calciumdxeventDetectAll');

bdeleteall = uicontrol('Style','pushbutton','Units','normalized','String','Delete events','Position',[.60 .05 .07 0.03],'FontSize',9,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; region.transients(1,num) = 1; hevPlotTrace;');

%popupmenu replacement for region.transient radio buttons, JBA 07/20/09
st = cell(1,5);
st{1} = 'none'; 
st{2} = 'sch';
st{3} = 'gdp';
st{4} = 'other';
st{5} = 'sink';
rtransients = uicontrol('Style','popupmenu','Units','normalized','String',st,'Position',[.70 .05 .07 0.03],'FontSize',9,...
    'BackgroundColor',[1 1 1]);

%Row 2 buttons
bdetect99 = uicontrol('Style','pushbutton','Units','normalized','String','Artifact rem w/ FFT notch','Position',[.10 .01 .07 0.03],'FontSize',9,...
    'Callback','calciumdxeventFFTnotchFilter');

%This following one is the good one
bdetect2 = uicontrol('Style','pushbutton','Units','normalized','String','Detect Artifacts','Position',[.20 .01 .07 0.03],'FontSize',9,...
    'Callback','calciumdxeventDetectArtifacts');

showArtifacts=0;
btoggle1 = uicontrol('Style','togglebutton','Units','normalized','String','Show Artifacts?','Position',[.30 .01 .07 0.03],'FontSize',9,...
    'Callback','hevArtifactToggle');

btdetect3=uicontrol('Style','pushbutton','Units','normalized','String','Manual Peak Find','Position',[.50 .01 .07 0.03],'FontSize',9,...
    'Callback','calciumdxManualPeaks');

showStimuli=0;
btStimImport=uicontrol('Style','togglebutton','Units','normalized','String','Show stimuli','Position',[.60 .01 .07 0.03],'FontSize',9,...
    'Callback','hevShowStimulusPeriods');

%{
%This following one is the good one
bdetect2 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Filt','Position',[.30 .01 .14 0.03],'FontSize',9,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP(trSign*region.traces,num,''no''); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hevPlotTrace;');
    
%------------- 
%The following is the good one, but for all traces
bdetect3 = uicontrol('Style','pushbutton','Units','normalized','String','Detect all Filt','Position',[.50 .01 .14 0.03],'FontSize',9,...
    'Callback','for c = 1:size(tr,1); [s d] = calciumdxEvent_DetSingTrHP(trSign*region.traces,c,''no''); set(progtx,''String'',[''Detecting '' num2str(c) '' of '' num2str(size(nt,1))]); spk(c,:) = 0; dec(c,:) = 0; spk(c,s) = 1; dec(c,d) = 1; if region.transients(1,c) == 1 && sum(spk(c,:)) > 0; region.transients(1,c) = 4; end; end; set(progtx,''String'',''''); hevPlotTrace;');

%---------------
 %This following one is the good one for high frequency signals
bdetect4 = uicontrol('Style','pushbutton','Units','normalized','String','Detect high freq','Position',[.70 .01 .14 0.03],'FontSize',9,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP_HF(trSign*region.traces,num,''no'',trNew); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hevPlotTrace;');
%}


% radio(1) = uicontrol('Style','radiobutton','Units','normalized','String','none','Position',[.60 .05 .07 0.04],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(radio(1),''value'',1); set(radio(2),''value'',0); set(radio(3),''value'',0); set(radio(4),''value'',0); hevTransients;');
%             
% radio(2) = uicontrol('Style','radiobutton','Units','normalized','String','sch','Position',[.67 .05 .07 0.04],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(radio(1),''value'',0); set(radio(2),''value'',1); set(radio(3),''value'',0); set(radio(4),''value'',0); hevTransients;');
% 
% radio(3) = uicontrol('Style','radiobutton','Units','normalized','String','gdp','Position',[.74 .05 .07 0.04],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(radio(1),''value'',0); set(radio(2),''value'',0); set(radio(3),''value'',1); set(radio(4),''value'',0); hevTransients;');
% 
% radio(4) = uicontrol('Style','radiobutton','Units','normalized','String','other','Position',[.81 .05 .07 0.04],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(radio(1),''value'',0); set(radio(2),''value'',0); set(radio(3),''value'',0); set(radio(4),''value'',1); hevTransients;');

bsave = uicontrol('Style','pushbutton','Units','normalized','String','Save','Position',[.90 .05 .07 0.04],'FontSize',12,...
    'Callback','hevSave');

num = 1;
hevPlotTrace
figure(fig)