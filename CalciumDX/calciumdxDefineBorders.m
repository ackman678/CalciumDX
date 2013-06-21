%Read resolution data
region.image = a;
region.spaceres = str2double(get(inptsr,'string'));
region.timeres = str2double(get(inpttr,'string'));
region.tifftype = tifftype;

if strcmp('OME',tifftype)
    if str2double(frameAveraging)<2
        frameaveraging = 'false';
    else
        frameaveraging = 'true';
    end
else
    frameaveraging = '';
end


prompt = {'Animal type (wt, C57BL6, GAD67.GFP, etc):','Age of animal (P7, adult, etc):','Type of experiment (img or img.ephys):','Calcium dye used:','Brain area (SC.L, V1.R, ctx, etc):','Enter field factor (f1, f2, etc):','Enter y or n if z-artifact present:','Enter z-depth below pia (in um):','Enter type of anesthetic:','Enter inhal anesthetic percentage:','Comments:','Enter whether frames were averaged:','Enter framerate (Hz):','Enter no. scanned lines per frame:','Enter scan line period:','Enter optical zoom:'};
dlg_title = 'Input experimental parameters';
num_lines = 1;
% def = {'wt','','img','CaGreenDextran','SC.R','','n','','isoflurane','1.0','',frameaveraging,framerate,linesperframe,scanlineperiod,opticalzoom};
def = {'wt','','img','CaGreenDextran','SC.R','','n','','isoflurane','1.0','',num2str(frameaveraging),num2str(framerate),num2str(linesperframe),num2str(scanlineperiod),num2str(opticalzoom)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

region.animaltype = answer{1};
region.age= answer{2};
region.exptype= answer{3};
region.dye= answer{4};
region.brainarea = answer{5};
region.field= answer{6};
region.zartifact= answer{7};
region.zdepth= answer{8};
region.anesthetic= answer{9};
region.anesthpercent= str2double(answer{10});
region.comments= answer{11};
region.frameaveraging= answer{12};
region.framerate= str2double(answer{13});
region.linesperframe= str2double(answer{14});
region.scanlineperiod= str2double(answer{15});
region.opticalzoom= str2double(answer{16});

delete(res_title);
delete(txlabsr);
delete(inptsr);
delete(txlabtr);
delete(inpttr);
delete(bnext);
set(bopenimage,'enable','off');

%{
maxdim=max(size(region.image));
if strcmp(region.brainarea,'SC.R')
    def_ans1 = 'anteriormedial point = (Y,X)';
    def_ans2 = [maxdim 0];
elseif strcmp(region.brainarea,'SC.L')
    def_ans1 = 'anteriormedial point = (Y,X)';
    def_ans2 = [0 0];
else
    def_ans1 = 'anteriormedial point = (Y,X)';
    def_ans2 = [0 0];
end

prompt = {'Provide description of orientation coordinate in image:','orientation coordinate in image [Y]:', 'orientation coordinate in image [X]:'};
dlg_title = 'Input experimental parameters';
num_lines = 1;
def = {def_ans1,num2str(def_ans2(1)),num2str(def_ans2(2))};
answer = inputdlg(prompt,dlg_title,num_lines,def);

region.orientation.description = answer{1};
region.orientation.value = [str2double(answer{2}) str2double(answer{3})];
%}

imgfig = figure; imagesc(region.image)
[x y butt] = ginput(1);
x = round(x);
y = round(y);
if x < 0
    x = 0;
elseif x > size(region.image,2)
    x = size(region.image,2);
end
if y < 0
    y = 0;
elseif y > size(region.image,1)
    y = size(region.image,1);
end

def_ans1 = 'anteriormedial point = (Y,X); (Y,X) = row, cols';
def_ans2 = [y x];
prompt = {'Provide description of orientation coordinate in image:','orientation coordinate in image [Y]:', 'orientation coordinate in image [X]:'};
dlg_title = 'Input experimental parameters';
num_lines = 1;
def = {def_ans1,num2str(def_ans2(1)),num2str(def_ans2(2))};
answer = inputdlg(prompt,dlg_title,num_lines,def);

region.orientation.description = answer{1};
region.orientation.value = [str2double(answer{2}) str2double(answer{3})];

close(imgfig)


%Contour functions
bord_title = uicontrol('Style','text','Units','normalized','String','Regions','Position',[.87 .755 .11 0.03],'FontSize',12,'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
bord_add = uicontrol('Style','pushbutton','Units','normalized','String','Add','Position',[.90 .595 .05 .03],'FontSize',9, ...
    'Enable','on','Callback','calciumdxAddBorder');
bord_delete = uicontrol('Style','pushbutton','Units','normalized','String','Delete','Position',[.90 .555 .05 .03],'FontSize',9, ...
    'Enable','off','Callback','calciumdxDeleteBorder');
bnext = uicontrol('Style','pushbutton','Units','normalized','String','Next >>','Position',[.63 .02 .05 .03],'FontSize',9, ...
    'Enable','on','Callback','calciumdxNameRegions');

%Initial info
bord = [];
bhand = [];

calciumdxDetermineRegions
axes(imgax)
