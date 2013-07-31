function region = myReadImageJROIregionAdd(region, overwriteCoords, scalingFactor)
%A function to read in ROIs from ImageJ and write them into a calciumdx formatted 'region' matlab data structure
%USE: region = myReadImageJROIregionAdd(region, 'false', 2)
%region-- set to 'true' or 'false' for whether you want to delete the current(if existing) coords
%overwriteCoords-- set to 'true' or 'false' for whether you want to delete the current(if existing) coords
%scalingFactor-- needed in case your ImageJ ROI coords were set on a half size working image, not the full resolution of the raw data
%After running this function and saving your 'region' data structure it is best to use the Import Contours function in calciumdx to import your new region.coords, especially if you want to use existing ROIs.  This function does not currently match your region.contours to your region.coords and assign them the region.location index. That is why running the Import Contours (contourarraysetup.m) is best
%James B. Ackman 2012-08-28
if nargin < 1 || isempty(region)
    try
        load calciumdxprefs
    end
    if exist('pathname','var')
        [filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open',pathname);
        if ~isstr(filename)
            delete(gcf)
            clear
        end
        if exist('calciumdxprefs.mat','file') == 2, save(calciumdxprefs,'pathname','filename','-append'); else, save(calciumdxprefs,'pathname','filename'); end
    else
        [filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open');
        if ~isstr(filename)
            delete(gcf)
            clear
        end
        if exist('calciumdxprefs.mat','file') == 2, save(calciumdxprefs,'pathname','filename','-append'); else, save(calciumdxprefs,'pathname','filename'); end
    end
    
    fnm = [pathname filename];
    load(fnm)
end

if nargin < 2 || isempty(overwriteCoords); overwriteCoords = 'true'; end
if nargin < 3 || isempty(scalingFactor); scalingFactor = 1; end

if strcmp(overwriteCoords,'false')
    count = length(region.coords); %current coord length, since we are going to add new coords. Use to initialize counter
    sROI = myReadImageJROI;  %function to read in ImageJ formatted ROIs from a file. Based on ReadImageJROI.m from matlabcentral
    for i = 1:length(sROI)
        count = count + 1;
        region.coords{count} = sROI{i}.mnCoordinates * scalingFactor;
        region.name{count} = sROI{i}.strName;
    end
    
else
    count = 1;
    region.coords = {};
    region.name = {};
    sROI = myReadImageJROI;
    for i = 1:length(sROI)
        region.coords{count} = sROI{i}.mnCoordinates * scalingFactor;
        region.name{count} = sROI{i}.strName;
        count = count + 1;
    end
    
end

%save([pathname filename],'region')
