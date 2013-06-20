%JBA, Monday, April 28, 2008 11:32 AM
clear
global bbright bcontrast

[filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open')
if ~isstr(filename)
    delete(gcf)
    clear
end
fnm1 = [pathname filename];
load(fnm1)
region1=region;

[filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open', pathname)
if ~isstr(filename)
    delete(gcf)
    clear
end
fnm2 = [pathname filename];
load(fnm2)
region2=region;

numOld=1;
opengl neverselect;
sz = get(0,'screensize');
fig = figure('Name','HippoEvents','NumberTitle','off','toolbar','figure','position',[1 0.15*sz(4) sz(3) 0.7*sz(4)],'doublebuffer','on');


uicontrol('Style','text','Units','normalized','String','Image','Position',[.87 .955 .11 0.03],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text','units','normalized','string','Brightness','position',[.87 .88 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bbright = uicontrol('Style','slider','Units','normalized','Position',[.87 .86 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','on','Callback','HippoContrast');
uicontrol('Style','text','units','normalized','string','Contrast','position',[.87 .83 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bcontrast = uicontrol('Style','slider','Units','normalized','Position',[.87 .81 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','on','Callback','HippoContrast');

a1=region1.image;
imgax1 = subplot('position',[0.1234375 0.703571 0.203125 0.282143]);
imagesc(a1);
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap gray
% set(bzoom,'enable','on');
set(bbright,'enable','on');
set(bcontrast,'enable','on');
%zoom on;

a2=region2.image;
imgax2 = subplot('position',[0.6234375 0.703571 0.203125 0.282143]);
imagesc(a2);
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap gray
% set(bzoom,'enable','on');
set(bbright,'enable','on');
set(bcontrast,'enable','on');
%zoom on;


[maxy maxx] = size(a1);


HippoContrast;
tr1 = region1.traces;
nt1 = [];
for c = 1:size(tr1,1)
    nt1(c,:) = dfoverf(tr1(c,:))*100;
end
spk1 = zeros(size(nt1));
dec1 = zeros(size(nt1));
for c = 1:size(spk1,1)
    spk1(c,region1.onsets{c}) = 1;
    dec1(c,region1.offsets{c}) = 1;
end
xlimits1 = [0 size(nt1,2)+1];
trax1 = subplot('position',[0.03 0.4 0.84 0.240357]);
box on
set(gcf,'KeyPressFcn','moveFroTrToTr1')
%zoom on;

tr2 = region2.traces;
nt2 = [];
for c = 1:size(tr2,1)
    nt2(c,:) = dfoverf(tr2(c,:))*100;
end
spk2 = zeros(size(nt2));
dec2 = zeros(size(nt2));
for c = 1:size(spk2,1)
    spk2(c,region2.onsets{c}) = 1;
    dec2(c,region2.offsets{c}) = 1;
end
xlimits2 = [0 size(nt2,2)+1];
trax2 = subplot('position',[0.03 0.1 0.84 0.240357]);
%zoom on;
box on


uicontrol('Style','text','Units','normalized','String','Cell #','Position',[0.87 0.5 0.05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
txcellnum1 = uicontrol('Style','edit','Units','normalized','String','1','Position',[.94 0.5 0.05 0.04],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','num1=dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);set(gcf,''KeyPressFcn'',''moveFroTrToTr1'')');
bgoto1 = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.89 .43 .05 0.04],'FontSize',12,...
    'Callback','num1=dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);set(gcf,''KeyPressFcn'',''moveFroTrToTr1'')');

uicontrol('Style','text','Units','normalized','String','Cell #','Position',[0.87 0.2 0.05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
txcellnum2 = uicontrol('Style','edit','Units','normalized','String','1','Position',[.94 0.2 0.05 0.04],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','num2=dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);set(gcf,''KeyPressFcn'',''moveFroTrToTr2'')');
bgoto2 = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.89 .13 .05 0.04],'FontSize',12,...
    'Callback','num2=dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);set(gcf,''KeyPressFcn'',''moveFroTrToTr2'')');


showAllcell1 = uicontrol('Style','pushbutton','Units','normalized','String','show contours','Position',[.35 .90 .10 0.04],'FontSize',12,...
    'Callback','allH1=showAllCont(region1,imgax1);');
delAllcell1 = uicontrol('Style','pushbutton','Units','normalized','String','delete contours','Position',[.35 .80 .10 0.04],'FontSize',12,...
    'Callback','delete(allH1)');
Selcell1 = uicontrol('Style','pushbutton','Units','normalized','String','select cell','Position',[.35 .70 .10 0.04],'FontSize',12,...
    'Callback','cellFound1=selCell(region1,imgax1);txtNum1');

showAllcell2 = uicontrol('Style','pushbutton','Units','normalized','String','show contours','Position',[.50 .90 .10 0.04],'FontSize',12,...
    'Callback','allH2=showAllCont(region2,imgax2);');
% showAllcell2 = uicontrol('Style','pushbutton','Units','normalized','String','show contours','Position',[.50 .90 .10 0.04],'FontSize',12,...
%     'Callback','allH2=showAllCont(region1,imgax2);');
delAllcell2 = uicontrol('Style','pushbutton','Units','normalized','String','delete contours','Position',[.50 .80 .10 0.04],'FontSize',12,...
    'Callback','delete(allH2)');
Selcell2 = uicontrol('Style','pushbutton','Units','normalized','String','select cell','Position',[.50 .70 .10 0.04],'FontSize',12,...
    'Callback','cellFound2=selCell(region2,imgax2);delete(txcellnum2);txtNum2');

num1=dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);set(gcf,'KeyPressFcn','moveFroTrToTr1')
num2=dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);set(gcf,'KeyPressFcn','moveFroTrToTr2')


% 
%         shaperad(1) = uicontrol('Style','radiobutton','Units','normalized','String','Circle','Position',[.87 .145 .05 0.025],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',1,'Callback','set(shaperad(1),''value'',1); set(shaperad(2),''value'',0);');
%         shaperad(2) = uicontrol('Style','radiobutton','Units','normalized','String','Custom','Position',[.925 .145 .055 0.025],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(shaperad(1),''value'',0); set(shaperad(2),''value'',1);');