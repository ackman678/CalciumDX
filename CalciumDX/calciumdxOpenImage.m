%removed OME-TIFF selector menu, 2013-06-19 17:00:35 JBA

try
    load calciumdxprefs
end
if exist('pathname','var')
    [filename, pathname] = uigetfile({'*.tif'; '*.lsm'}, 'Choose image to open',pathname);
    if ~ischar(filename)
        return
    end
else
    [filename, pathname] = uigetfile({'*.tif'; '*.lsm'}, 'Choose image to open');
    if ~ischar(filename)
        return
    end
end
fnm = [pathname filename];
save(calciumdxprefs,'pathname', 'filename')

frameAveraging = '';
framerate = '';
linesperframe = '';
opticalzoom = '';
scanlineperiod = '';

% Now we use loci_tools.jar (make sure is present in local path) to use the
% Open Microscopy Environment (OME) java plugin to flexibly open the time
% series (should work with a variety microscope manufacturer formats.
% bfopen.m might have to be edited to ensure compatibility.

data = bfopen(fnm)
tifftype = 'OME';


% chanlabs = vertcat(data{1,1}{:,3});  %jba 2011-02-09
% numchan = max(chanlabs);   %jba 2011-02-09

numchan = str2double(data{1,2}.get('NumberOfChannels'));

%must fix first part of this if statement to use with the channel separator implemented in new bfopen.m and loci_tools.jar from 2011-06-07
%Fixed this script and bfopen on 2011-02-09 to handle timeseries with more than one channel. Must select the one channel to analyse for now.
if numchan > 1  %jba 2011-02-09
    numframes = length(find(chanlabs < 2));
    sz = size(data{1, 1}{1,1});
    prompt = {['File has ' num2str(numchan) ' channels. Enter channel number']};
    dlg_title = 'Enter chan no.';
    num_lines = 1;
    def = {'2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    chanIdx = str2double(answer{1});
    series1=zeros([sz numframes],'uint16');
    IND = find(chanlabs == chanIdx);
    for i=1:length(IND)
        series1(:,:,i)=data{1, 1}{IND(i),1};
    end
else
    numframes = size(data{1, 1},1);
    sz = size(data{1, 1}{1,1});
    series1=zeros([sz numframes],'uint16');   %jba 2011-06-07. Directlly specified'uint16' class
    for i=1:numframes
        series1(:,:,i)=data{1, 1}{i,1};
    end
end
%jba 2011-06-07.  added following lines to clear out large java image series array.
tmp = data{1,2};
clear data;
data{1,1} = {};
data{1,2} = tmp;

umpx = data{1,2}.get('micronsPerPixel_XAxis'); %umperpx, but gives double the state value instead half like it should...
frameAveraging = data{1,2}.get('frameAveraging');
framerate = data{1,2}.get('framerate'); %Hz
frameperiod = str2double(data{1,2}.get('framePeriod'));

%the following code worked well for differentiating whole frameAvergaing for single channel data, but not for multichannel data for some reason. So commented out on 2011-02-09
%{
    if strcmp(frameAveraging,'0')
        frameAveraging = '1';
        frameperiod = str2double(data{1,2}.get('framePeriod')) * str2double(frameAveraging);
    else
        frameperiod = str2double(data{1,2}.get('framePeriod')) * str2double(frameAveraging); %sec
    end
%}
linesperframe = data{1,2}.get('fullImageLinesPerFrame');
scanlineperiod = data{1,2}.get('scanlinePeriod');
opticalzoom = data{1,2}.get('opticalZoom');

umpxReal = str2double(umpx)/str2double(opticalzoom);  %here we calculate the true spatial resolution based on the opticalzoom parameter


set(inptsr,'String',umpxReal)
set(inpttr,'String',frameperiod)


% pr = zeros(1,numframes,3);
region = [];

if numframes == 1
    a = double(imread(fnm));
    button = questdlg({'The file contains one frame.','How was it obtained?'},'Frame information','Average','Maximum','First frame','Average');
    if strcmp(button,'Average')
        region.frametype = 'average';
    elseif strcmp(button,'Maximum')
        region.frametype = 'maximum';
    elseif strcmp(button,'First frame')
        region.frametype = 'first';
    end
else
    button = questdlg({['The file contains ' num2str(numframes) ' frames.'],'Choose the operation to perform.'},'Frame information','Average','Maximum','First frame','Average');
    if strcmp(button,'Average')
        region.frametype = 'average';
        
        
        %Should implement imshow(mat2gray(A)) here so user can visually pick some subsequent frames that have no movement for avg image
        %----------------------------
        %--decide whether to average all frames (default), or a subset
        %(good for movies that shake a little bit in the xy direction
        %to outline clear cells
        prompt = {'Start','Last'};
        dlg_title = 'Average frames';
        num_lines = 1;
        def = {'1',num2str(numframes)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        numframesAvg = str2double(answer{1});
        baseframes = [str2double(answer{1}) str2double(answer{2})];
        %----------------------------
        
        %old way avg method
        %             for c = 1:numframes
        %                 %             if mod(c,20)==1
        %                 %                 prg = subplot('position',[.87 .575 .11 0.025]);
        %                 %                 imagesc(pr);
        %                 %                 axis off
        %                 %                 drawnow;
        %                 %                 delete(prg);
        %                 %             end
        %                 h = waitbar(c/numframes);
        %                 a = a + series1{c,1};
        %                 %             pr(1,c,1) = 1;
        %             end
        
        %new avg method
        a=mean(series1(:,:,baseframes(1):baseframes(2)),3);
        %             close(h)
        %             a = a / numframes;
    elseif strcmp(button,'Maximum')
        region.frametype = 'maximum';
        %old method
        %             a = zeros(size(series1{1,1}));
        %             for c = 1:numframes
        %                 %             if mod(c,20)==1
        %                 %                 prg = subplot('position',[.87 .575 .11 0.025]);
        %                 %                 imagesc(pr);
        %                 %                 axis off
        %                 %                 drawnow;
        %                 %                 delete(prg);
        %                 %             end
        %                 h = waitbar(c/numframes);
        %                 a(:,:,2) = series1{c,1};
        %                 a(:,:,1) = max(a,[],3);
        %                 %             pr(1,c,1) = 1;
        %             end
        %             close(h)
        %             a = double(a(:,:,1));
        
        a=max(series1,3);
    elseif strcmp(button,'First frame')
        region.frametype = 'first';
        %             a = series1{c,1};  %old method
        a = series1(:,:,1);
    end
end

imgax = subplot('position',[0.02 0.02 0.82 0.96]);
imagesc(a);
hold on

set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on

colormap gray

set(bzoom,'enable','on');
set(bbright,'enable','on');
set(bcontrast,'enable','on');
set(bnext,'enable','on');
set(inptsr,'enable','on');
set(inpttr,'enable','on');


[maxy maxx] = size(a);

zoom on;
calciumdxContrast;