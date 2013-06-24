%Read resolution data
region.image = a;
region.spaceres = str2double(get(inptsr,'string'));
region.timeres = str2double(get(inpttr,'string'));
region.tifftype = tifftype;

if strcmp('OME',tifftype) && ~isempty(frameAveraging)
    if str2double(frameAveraging)<2
        frameaveraging = 'false';
    else
        frameaveraging = 'true';
    end
else
    frameaveraging = '';
end

region.animaltype = '';
region.age= '';
region.exptype= '';
region.dye= '';
region.brainarea = '';
region.field= '';
region.zartifact= '';
region.zdepth= '';
region.anesthetic= '';
region.anesthpercent= '';
region.comments= '';

region.frameaveraging= frameaveraging;
region.framerate= str2double(framerate);
region.linesperframe= str2double(linesperframe);
region.scanlineperiod= str2double(scanlineperiod);
region.opticalzoom= str2double(opticalzoom);


btProps1 = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Exp props','Position',[.87 .81-0.03 .04 0.025],'FontSize',9,'Callback','dxInputParamsExp');
btProps2 = uicontrol(fig,'Style','pushbutton','Units','normalized','String','2P props','Position',[.92 .81-0.03 .04 0.025],'FontSize',9,'Callback','dxInputParams2P');

delete(res_title);
delete(txlabsr);
delete(inptsr);
delete(txlabtr);
delete(inpttr);
delete(bnext);
set(bopenimage,'enable','off');

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
bord_title = uicontrol('Style','text','Units','normalized','String','Regions','Position',[.87 .735 .11 0.03],'FontSize',12,'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
bord_add = uicontrol('Style','pushbutton','Units','normalized','String','Add','Position',[.90 .595 .05 .03],'FontSize',9, ...
    'Enable','on','Callback','calciumdxAddBorder');
bord_delete = uicontrol('Style','pushbutton','Units','normalized','String','Delete','Position',[.90 .555 .05 .03],'FontSize',9, ...
    'Enable','off','Callback','calciumdxDeleteBorder');
bnext = uicontrol('Style','pushbutton','Units','normalized','String','Next >>','Position',[.87 .05 .05 .03],'FontSize',9, ...
    'Enable','on','Callback','calciumdxNameRegions');
    
%Initial info
bord = [];
bhand = [];

calciumdxDetermineRegions
axes(imgax)
