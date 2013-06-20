function getWaveMasks(filelist,region,wavesDBfnm,series1)
%getWaveMasks({filename},region,[],series1);  %for a single file
%getWaveMasks(filelist);  %for a list of files
%(1) table with desired filenames (space delimited txt file, with full filenames for region.mat files and tiff movies in first and second columns respectively)
%files.txt should have filenames in first column.
%can have extra columns with descriptor/factor information for the file. This will be the rowinfo that is attached to each measure observation in the following script.
%filelist = readtext('files.txt',' '); %grab readtext.m file script from matlab central
%or
%(2) a single filename (filename of your region .mat file) as a cell array, i.e.  {filename}
%James B. Ackman 2011-11-28

%LOAD FILE
if nargin< 3 || isempty(wavesDBfnm); 
wavesDBfnm = 'wavesDB.mat'; 
%wavesDB = []; 
filename = []; tiffFilename = []; orientation = []; regionname = []; regioncoords = []; wavenumber = []; maskparams = []; imagedata = []; 
%save(wavesDBfnm,'wavesDB','-v7.3'); 

if exist(fullfile(cd, wavesDBfnm),'file') == 0
save(wavesDBfnm,'filename', 'tiffFilename', 'orientation', 'regionname', 'regioncoords', 'wavenumber', 'maskparams', 'imagedata', '-v7.3');
else
error([wavesDBfnm ' already exists!'])
end
end

if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end
%results={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency'};
fnms=filelist(:,1);
if nargin < 4 || isempty(series1); tiffFnms=filelist(:,2); end

h = waitbar(0,'Please wait...');
for j=1:numel(fnms)
    waitbar(j/numel(fnms),h);
    
    if loadfile > 0
        myMATfile=load(fnms{j});
        region=myMATfile.region;
    end
    [pathstr, name, ext] = fileparts(fnms{j});
    rowinfo = {[name ext]};  %2011-07-11 jba
    %     rowinfo = filelist(j,:);
    sprintf(fnms{j})
    filename = fnms{j};
    %open tiff array
    if nargin < 4 || isempty(series1)
        tiffFilename = tiffFnms{j};
        [data, series1] = myOpenOMEtiff(tiffFilename);
    else
        tiffFilename = fnms{j};  %TESTING
    end
    
%    if isempty(wavesDB)
%        filecounter = 0;
%    else
%        filecounter = length(wavesDB);
%    end
    %-----------------------------------------------------------------------------------------
    %get masks by region
    locationMarkers = unique(region.location);
    for i = 1:length(locationMarkers)
        locationIndex = locationMarkers(i);
		wavemasks = getWaveMasksByLocation(region,locationIndex,series1,filename,tiffFilename);
%		if isempty(wavesDB)
%			filecounter = 0;
%		else
%			filecounter = length(wavesDB);
%		end
		
		
        if isempty(wavemasks)
            continue
        else
%			filecounter = filecounter +1;			
%            wavesDB(filecounter+1).filename = tiffFilename;
%            wavesDB(filecounter+1).orientation = region.orientation;
%            wavesDB(filecounter+1).regionname(i).waves = wavemasks;
%            wavesDB(filecounter+1).regionname(i).name = region.wavedata{locationIndex}.name;
%            wavesDB(filecounter).waves.filename = tiffFilename;
%            wavesDB(filecounter).waves.orientation = region.orientation;
%			 wavesDB(filecounter).waves.region = region.wavedata{locationIndex}.name;            
%            wavesDB(filecounter).waves.data = wavemasks;
            
%            wavesDB = [wavesDB wavemasks];
			matObj = matfile(wavesDBfnm,'Writable',true);   %open wavesDB for writing to file. wavesDBfnm must be full file name or relative to working dir.
%			szDB = size(matObj,'wavesDB',2);
			szDB = size(matObj,'filename',2);
%			names = matObj.wavesDB(1).filename
			sz_wavemasks = numel(wavemasks.filename);
			
			if szDB == 0;
%				matObj.wavesDB = wavemasks;   %append new waves to existing wavesDB
				matObj.filename = wavemasks.filename;   %append new waves to existing wavesDB
				matObj.tiffFilename = wavemasks.tiffFilename;   %append new waves to existing wavesDB				
				matObj.orientation = wavemasks.orientation;   %append new waves to existing wavesDB
				matObj.regionname = wavemasks.regionname;   %append new waves to existing wavesDB
				matObj.regioncoords = wavemasks.regioncoords;   %append new waves to existing wavesDB
				matObj.wavenumber = wavemasks.wavenumber;   %append new waves to existing wavesDB
				matObj.maskparams = wavemasks.maskparams;   %append new waves to existing wavesDB
				matObj.imagedata = wavemasks.imagedata;   %append new waves to existing wavesDB
			else
%				matObj.wavesDB(:,szDB+1:szDB+numel(wavemasks)) = wavemasks;   %append new waves to existing wavesDB
				matObj.filename(1,szDB+1:szDB+sz_wavemasks) = wavemasks.filename;   %append new waves to existing wavesDB
				matObj.tiffFilename(1,szDB+1:szDB+sz_wavemasks) = wavemasks.tiffFilename;   %append new waves to existing wavesDB				
				matObj.orientation(1,szDB+1:szDB+sz_wavemasks) = wavemasks.orientation;   %append new waves to existing wavesDB
				matObj.regionname(1,szDB+1:szDB+sz_wavemasks) = wavemasks.regionname;   %append new waves to existing wavesDB
				matObj.regioncoords(1,szDB+1:szDB+sz_wavemasks) = wavemasks.regioncoords;   %append new waves to existing wavesDB
				matObj.wavenumber(1,szDB+1:szDB+sz_wavemasks) = wavemasks.wavenumber;   %append new waves to existing wavesDB
				matObj.maskparams(1,szDB+1:szDB+sz_wavemasks) = wavemasks.maskparams;   %append new waves to existing wavesDB
				matObj.imagedata(1,szDB+1:szDB+sz_wavemasks) = wavemasks.imagedata;   %append new waves to existing wavesDB
			end
        end
    end
    
    %-----------------------------------------------------------------------------------------
    %data=results;
end
close(h)
end




function wavemasks = getWaveMasksByLocation(region,locationIndex,series1,filename,tiffFilename)
if isfield(region.wavedata{locationIndex},'waveonsets')
    if ~isempty(region.wavedata{locationIndex}.waveonsets) && (strcmp(region.wavedata{locationIndex}.name,'SC.R') | strcmp(region.wavedata{locationIndex}.name,'SC.L'))
        nwaves = numel(region.wavedata{locationIndex}.waveonsets);
        locationMarkers = unique(region.location);
        sz = size(series1);
		regionMask1 = poly2mask(region.coords{locationIndex}(:,1),region.coords{locationIndex}(:,2),sz(1),sz(2));
		coordXYvalues = vertcat(region.coords{[locationMarkers]});
		minXcoord = unique(min(coordXYvalues(:,1)));
		maxXcoord = unique(max(coordXYvalues(:,1)));
		
        hbar = waitbar(0,'Please wait...');
        for j = 1:nwaves
        	clear A B BWarray
            waitbar(j/nwaves,hbar);
            
            waveONidx=region.wavedata{locationIndex}.waveonsets(j);
            waveOFFidx=region.wavedata{locationIndex}.waveoffsets(j);
            A = double(series1(:,:,waveONidx:waveOFFidx));  %make raw image array
            sz = size(A);  %array size
            
            for fr = 1:sz(3)
                A(:,:,fr) = (A(:,:,fr) - region.image)./region.image;   %replace each raw image with the normalized dF/F representation
            end
            
            %----Make Binary Mask Movie of the waveframes--------------
%             figure;
            %clear M
            
            %	figure; imshow(regionMask1);
            cropSZ=region.orientation.value(1);
            gaussSmoothSigma = 6;
            graythreshlevel = [];
            absoluteMinimumDfoverfValue = [];
            for fr = 1:sz(3)
                img1 = gaussSmooth(A(:,:,fr),gaussSmoothSigma,'same');  %accentuate big blobby objects at 6 sigma, about 50-100um diam.
                
 %               fr1=img1(1:cropSZ,:);   %use these lines for cropping and flipping the masks for comparisons between hemispheres.
%                fr2=img1(cropSZ+1:sz(1),:);
%                fr2=flipud(fr2);
                
                level = graythresh(img1);
                BW = im2bw(img1,level);
                fr1ind = intersect(find(regionMask1),find(BW));
%				BW2 = zeros(size(BW),'uint8'); %to save the binary mask
%				BW2(fr1ind) = BW(fr1ind);
%				BWarray(:,:,fr) = BW2(:,minXcoord:maxXcoord);

				BW2 = zeros(size(BW));	%to save the mask with the gaussSmooth intensity values instead
				BW2(fr1ind) = img1(fr1ind);
				tmp = BW2(:,minXcoord:maxXcoord);  
				tmp2= (tmp + abs(min(tmp(:)))).*2^16;  %scale data to fit within uint16 range [0,2^16)				
				BWarray(:,:,fr) = uint16(tmp2);
                
                %	imshow(BW2,[])  %TESTING
                %	M(fr) = getframe;
                

%				tmp = A(:,minXcoord:maxXcoord,fr);  %to save the raw dfoverf image
				tmp = img1(:,minXcoord:maxXcoord);  %to save the GaussSmooth image
				minValue = abs(min(tmp(:)));
				tmp2= (tmp + minValue).*2^16;  %scale data to fit within uint16 range [0,2^16)
                B(:,:,fr) = uint16(tmp2);    %store as uint16
                graythreshlevel = [graythreshlevel; level];
                absoluteMinimumDfoverfValue = [absoluteMinimumDfoverfValue; minValue];
            end

            %-------add waves to db-----------
%			wavemasks.dfoverf = B;
%            wavemasks.binarymask = BWarray;
            newcoords = region.coords{locationIndex};   %to shift the region.coord polygon coords outline we will save since the image is now cropped
			newcoords(:,1) = newcoords(:,1) - (minXcoord-1);
            
            wavemasks.filename(j).filename = filename;
            wavemasks.tiffFilename(j).tiffFilename = tiffFilename;            
            wavemasks.orientation(j).orientation = region.orientation;
			wavemasks.regionname(j).regionname = region.wavedata{locationIndex}.name;
			wavemasks.regioncoords(j).regioncoords = newcoords;            
			wavemasks.wavenumber(j).wavenumber = j;
			wavemasks.maskparams(j).gaussSmoothSigma = gaussSmoothSigma;
			wavemasks.maskparams(j).graythreshlevel = graythreshlevel;
			wavemasks.maskparams(j).absoluteMinimumDfoverfValue = absoluteMinimumDfoverfValue;
            wavemasks.imagedata(j).dfoverf = B;
            wavemasks.imagedata(j).wavemask = BWarray;            
            
%            wavemasks(j).dfoverf = B;
%            wavemasks(j).binarymask = BWarray;
            %clear A B BWarray
            %movie(M,30)
            %vidObj = VideoWriter(['xcorrn_' datestr(now,'yyyymmdd-HHMMSS') '.avi'])
            %open(vidObj)
            %for i =1:numel(M)
            %writeVideo(vidObj,M(i))
            %end
            %close(vidObj)
        end
        close(hbar)
        
    else
        wavemasks = [];
    end
end
end
