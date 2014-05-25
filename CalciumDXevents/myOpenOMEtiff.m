function [data, series1, filename] = myOpenOMEtiff(fnm)
%Open time series tiff stacks
%uses the bfopen OME tiff matlab function from openmicroscopy.org
%bfopen uses the loci_tools.jar OME library to open all types of microscopy data
%returns any metadata and the image stack as a uint16 array
%loci_tools.jar and bfopen.m has to be in your matlab search path
%
% USAGE
%  myOpenOMEtiff;
%
% INPUTS
%  fnm           - optional full filename
%
% OUTPUTS
%  data           - structure containing metadata for the image series
%  series1        - image series
%
%James Ackman, 2011-11-10

if nargin < 1 || isempty(fnm)
    
    if exist('calciumdxprefs.mat','file') == 2
        load('calciumdxprefs')
    else
        pathname = pwd;
    end
    
    if exist('pathname','var')
        [filename, pathname] = uigetfile({'*.tif'}, 'Choose image to open',pathname);
        if ~ischar(filename)
            return
        end
    else
        [filename, pathname] = uigetfile({'*.tif'}, 'Choose image to open');
        if ~ischar(filename)
            return
        end
    end
    fnm = [pathname filename];
    save('calciumdxprefs.mat', 'pathname','filename')
end

%     tic
    data = bfopen(fnm)
%     toc

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

% assignin('base', 'data',data)
% assignin('base', 'series1',series1)