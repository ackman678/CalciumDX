%calciumdx.m
%for analysis of calcium imaging movies (either 2P or CCD based)
%2013-06-24 miscellaneous bug fixes
%2013-06-21 miscellaneous bug fixes
%2013-06-20 10:35:50 started fork and reorganization by James Ackman for upload to github as shared project editing. Many TODOs are needed including reworking as object oriented.
% JBA edited from Nov 2011- Mar 2012 for importing and aligning contours and coords, and indexing/ contour out of image bounds error
%Added bug fixes at v4.1 by James Ackman from Nov 2011- Mar 2012 for importing and aligning contours and coords, and indexing/ contour out of image bounds error
%reworked by James Ackman at v4.0 between October 2010 - October 2011 in Crair lab at Yale as hippocalciumdextran for analysis of retinal wave movies (Ackman et al., Nature 2012)
%clear all contours, import contours, align all contours, align single contours added by James Ackman July 19, 2009
%capability for opening open microscopy environment tiff formats (http://www.loci.wisc.edu/ome/ome-tiff.html) added by James Ackman July 17, 2009.
%repeat delete added to calciumdxManualDelete.m by James Ackman 17/01/2007
%Default spatial and temporal values edited 17/01/2007
%2006-2007 forked from an earlier project version named hippo first created by D.Aronov and R. Cossart at INMED.

calciumdxFullfilename = mfilename('fullpath');
[calciumdxpath, ~, ~] = fileparts(calciumdxFullfilename);
matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxbackupLocation = fullfile(matlabUserPath,'calciumdxbackup.mat');
calciumdxprefs = fullfile(matlabUserPath,'calciumdxprefs.mat');


save(calciumdxbackupLocation, 'matlabUserPath')
%clear;
cd(calciumdxpath)


versionnum = '4.1.2';

if isdir('TraceReaders')
    
    %opengl neverselect;
	%fig = figure('Name',['calciumdx ' versionnum],'NumberTitle','off','MenuBar','none','doublebuffer','on','units','normalized','CloseRequestFcn','calciumdxFigClose','position',[0 .08/3 1 2.86/3]);
    sz = get(0,'screensize');
    fig = figure('Name','calciumdx ','NumberTitle','off','MenuBar','none','CloseRequestFcn','calciumdxFigClose','position',[1 sz(4) sz(3) sz(4)],'doublebuffer','on');
	old_units = get(fig,'Units');
	set(fig,'Units','normalized');
    
    %Image functions
    uicontrol('Style','text','Units','normalized','String','Image','Position',[.87 .955 .11 0.03],'FontSize',12,'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
    bopenimage = uicontrol('Style','pushbutton','Units','normalized','String','Open','Position',[.87 .91 .05 .03],'FontSize',9, ...
        'Callback','calciumdxOpenImage');
    bzoom = uicontrol('Style','pushbutton','Units','normalized','String','Zoom','Position',[.93 .91 .05 .03],'FontSize',9, ...
        'Callback','zoom on','Enable','off');
    uicontrol('Style','text','units','normalized','string','Brightness','position',[.87 .88 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
    bbright = uicontrol('Style','slider','Units','normalized','Position',[.87 .86 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
        'Enable','off','Callback','calciumdxContrast');
    uicontrol('Style','text','units','normalized','string','Contrast','position',[.87 .83 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
    bcontrast = uicontrol('Style','slider','Units','normalized','Position',[.87 .81 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
        'Enable','off','Callback','calciumdxContrast');
    
    res_title = uicontrol('Style','text','Units','normalized','String','Resolution','Position',[.87 .755 .11 0.03],'FontSize',12,'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
    txlabsr = uicontrol('Style','text','Units','normalized','String','Spatial (µm/pixel)','Position',[.87 .715 .11 0.02],'FontSize',9,...
        'BackgroundColor',[.8 .8 .8],'HorizontalAlignment','left');
    inptsr = uicontrol('Style','edit','Units','normalized','String','','Position',[.87 .715-0.0275 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left','enable','off');
    txlabtr = uicontrol('Style','text','Units','normalized','String','Temporal (sec/frame)','Position',[.87 .715-0.0275-0.025 .11 0.02],'FontSize',9,...
        'BackgroundColor',[.8 .8 .8],'HorizontalAlignment','left');
    inpttr = uicontrol('Style','edit','Units','normalized','String','','Position',[.87 .715-2*0.0275-0.025 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left','enable','off');
    
    bnext = uicontrol('Style','pushbutton','Units','normalized','String','Next >>','Position',[.87 .05 .05 .03],'FontSize',9, ...
        'Enable','off','Callback','calciumdxDefineBorders');
    
else
    errordlg('Please change working directory to the folder containing calciumdx.m','Program Error');
end
