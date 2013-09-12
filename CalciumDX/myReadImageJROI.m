function output = myReadImageJROI
%A function to read in ROIs from ImageJ, where these ImageJ roi file is typically a '.zip' archive for a set of ROIs or a single ROI with the extension '.roi'
%Uses ReadImageJROI.m from matlabcentral
%James B. Ackman 2012-08-28
[filename, pathname] = uigetfile({'*.zip' '*.roi'}, 'Choose data file to open');
% tempdir=pwd;
% cd ReadImageJROI;
[sROI] = ReadImageJROI([pathname filename]);
% cd(tempdir);
output = sROI;