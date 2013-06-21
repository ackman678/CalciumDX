if strcmp(get(det_view,'enable'),'off')
    for c = 1:length(handl)
        delete(handl{c});
    end
    cn = cell(1,length(region.name));
    centr = cell(1,length(region.name));
    areas = cell(1,length(region.name));
    handl = cell(1,length(region.name));
    handlCoord = cell(1,length(region.name));
    set(det_view,'enable','on');
    
    thres = 30*ones(1,length(region.name)); %previously default of 15
    old_thres = inf*ones(1,length(region.name));
    lowar = 5*ones(1,length(region.name)); %previously default of 25
    highar = repmat(inf,1,length(region.name));
    pilim = 4*ones(1,length(region.name));
    isadjust = zeros(1,length(region.name));
    isdetected = zeros(1,length(region.name));
    ishid = 1;
    isimported = 0;
    
    if islocal == 0
        
        %---START new 07/19/09       
        %This button will clear your current contours.
        btclear = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Clear all','Position',[.02 .98 .1 0.03],'FontSize',9,...
            'Callback','cn{num} = []; calciumdxDrawCells;');
        %This button will import contour data from an existing .mat file.
        btimport = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Import Contours','Position',[.02 .94 .1 0.03],'FontSize',9,...
            'Callback','[tfilename, tpathname] = uigetfile(); tfnm=fullfile(tpathname,tfilename); tmp = load(tfnm); region.contours=tmp.region.contours; region.name=tmp.region.name; region.coords=tmp.region.coords; region.location=tmp.region.location; region.cutoff=tmp.region.cutoff; region.lowarea=tmp.region.lowarea; region.higharea=tmp.region.higharea; region.isdetected=tmp.region.isdetected; region.pilimit=tmp.region.pilimit; region.isadjusted=tmp.region.isadjusted; isimported = 1; cn = cell(1,length(region.name)); centr = cell(1,length(region.name)); areas = cell(1,length(region.name)); handl = cell(1,length(region.name)); handlCoord = cell(1,length(region.name)); contourarraysetup; lowar = zeros(1,length(region.name)); highar = repmat(inf,1,length(region.name)); thres = 30*ones(1,length(region.name)); pilim = 4*ones(1,length(region.name)); cl = hsv(length(region.name)); calciumdxDrawCells;');
            
        %This button will register the contours over the cells. User clicks left-mouse button on desired location then clicks on original location. Will automatically redraw and assign contours to the newly specified location.
        btalign = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Align All','Position',[.02 .90 .1 0.03],'FontSize',9,...
            'Callback','[x,y] = ginput(2); dx = x(1) - x(2); dy = y(1) - y(2); for c = 1:length(cn{num}); cn{num}{c}(:,1) = cn{num}{c}(:,1) - dx; cn{num}{c}(:,2) = cn{num}{c}(:,2) - dy; end; calciumdxDrawCells;');     
       
        btalign2 = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Align Single','Position',[.02 .86 .1 0.03],'FontSize',9,...
            'Callback','[x,y] = ginput(2); matr = []; for c = 1:length(cn); for d = 1:length(cn{c}); if polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) > lowar(c) & polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) < highar(c); matr = [matr; [calciumdxCentroid(cn{c}{d}) c d]]; end; end; end; dst = sum((matr(:,1:2)-repmat([x(1) y(1)],size(matr,1),1)).^2,2); [dummy i] = min(dst); dx = x(1) - x(2); dy = y(1) - y(2); cn{num}{i}(:,1) = cn{num}{i}(:,1) - dx; cn{num}{i}(:,2) = cn{num}{i}(:,2) - dy; calciumdxDrawCells;');
        
        btalign3 = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Align All + Coords','Position',[.02 .82 .1 0.03],'FontSize',9,'Callback','calciumdxRealignCoordsContours');
        %---END new 07/19/09
        
        txlab = uicontrol(fig,'Style','text','Units','normalized','String',region.name{num},'Position',[.87 .49 .11 0.025],'FontSize',10,'FontWeight','Bold',...
            'BackgroundColor',cl(num,:));
        dummyp(1) = uicontrol(fig,'Style','text','Units','normalized','String','Cutoff','Position',[.87 .4625 .11 0.02],'FontSize',9,...
            'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
        txthres = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(thres(num)),'Position',[.93 .46 .05 0.025],'FontSize',9,...
            'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
        dummyp(2) = uicontrol(fig,'Style','text','Units','normalized','String','Min area','Position',[.87 .4325 .11 0.02],'FontSize',9,...
            'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
        txarlow = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(lowar(num)),'Position',[.93 .43 .05 0.025],'FontSize',9,...
            'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
        dummyp(3) = uicontrol(fig,'Style','text','Units','normalized','String','Max area','Position',[.87 .4025 .11 0.02],'FontSize',9,...
            'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
        txarhigh = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(highar(num)),'Position',[.93 .40 .05 0.025],'FontSize',9,...
            'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
        
        cmnd = ['thres(num) = str2num(get(txthres,''string'')); lowar(num) = str2num(get(txarlow,''string'')); highar(num) = str2num(get(txarhigh,''string'')); pilim(num) = str2num(get(txpilim,''string''));'];
        btdetect = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Detect!','Position',[.87 .355 .05 0.03],'FontSize',9,...
            'Callback',[cmnd 'calciumdxFindCells']);
        bthide = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Hide','Position',[.93 .355 .05 0.03],'FontSize',9,...
            'Callback','ishid=1-ishid; calciumdxHide');
        
        dummyp(4) = uicontrol(fig,'Style','text','Units','normalized','String','Pi limit','Position',[.87 .3075 .11 0.02],'FontSize',9,...
            'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
        txpilim = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(pilim(num)),'Position',[.93 .305 .05 0.025],'FontSize',9,...
            'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
        btfindbad = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Find','Position',[.87 .26 .05 0.03],'FontSize',9,...
            'Callback','calciumdxFindBad');
        btadjust = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Adjust','Position',[.93 .26 .05 0.03],'FontSize',9,...
            'Callback','calciumdxAdjust');
        
        btprev = uicontrol(fig,'Style','pushbutton','Units','normalized','String','<< Prev','Position',[.87 .205 .05 0.03],'FontSize',9,...
            'Callback',[cmnd 'ishid=1; calciumdxHide; num=mod(num+length(region.name)-2,length(region.name))+1; calciumdxInputParams;']);
        btnext = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .205 .05 0.03],'FontSize',9,...
            'Callback',[cmnd 'ishid=1; calciumdxHide; thres(num) = str2num(get(txthres,''string'')); num=mod(num,length(region.name))+1; calciumdxInputParams;']);
        
        dummyp(5) = uicontrol(fig,'Style','text','Units','normalized','String','Manual add ROIs','Position',[.87 .1725 .11 0.02],'FontSize',9,...
            'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
%         shaperad(1) = uicontrol(fig,'Style','radiobutton','Units','normalized','String','Add Region','Position',[.87 .145 .05 0.025],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',1,'Callback','set(shaperad(1),''value'',1); set(shaperad(2),''value'',0);');
%         shaperad(2) = uicontrol(fig,'Style','radiobutton','Units','normalized','String','Add','Position',[.87 .120 .055 0.025],'FontSize',9,...
%             'BackgroundColor',[.8 .8 .8],'Value',0,'Callback','set(shaperad(1),''value'',0); set(shaperad(2),''value'',1);'); 

        allowRegionOverlap=0;
        btoggle1 = uicontrol('Style','togglebutton','Units','normalized','String','Allow Region overlap?','Position',[.87 .145 .05 0.03],'FontSize',9,...
            'Callback','hevAllowOverlapToggle');

        btadd = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Add','Position',[.87 .08 .05 0.03],'FontSize',9,...
            'Callback','calciumdxManualAdd');        
        btdelete = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Delete','Position',[.93 .08 .05 0.03],'FontSize',9,...
            'Callback','calciumdxManualDelete','enable','off');
        
        btnextscr = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 0.03],'FontSize',9,...
            'Callback','calciumdxReadTraceParam','enable','on');
        islocal = 1;
    end    
end

set(txlab,'String',region.name{num},'BackgroundColor',cl(num,:));
set(txthres,'String',num2str(thres(num)));
set(txarlow,'String',num2str(lowar(num)));
set(txarhigh,'String',num2str(highar(num)));
set(txpilim,'String',num2str(pilim(num)));

zoom on

% btcalciumdxsavebackup = uicontrol('Style','pushbutton','Units','normalized','String','Save backup','Position',[.02 .04 .05 .03],'FontSize',9, 'Callback','calciumdxSaveBackup');
% btcalciumdxloadbackup = uicontrol('Style','pushbutton','Units','normalized','String','Load backup','Position',[.02 .01 .05 .03],'FontSize',9, 'Callback','calciumdxLoadBackup');

% hBar = waitbar(0,'Please wait...');
% currentProgram = 'calciumdxInputParams.m';
%save('tmpcalciumdxBackup','region','cn')
% delete(hBar)