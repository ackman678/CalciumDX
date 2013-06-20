function [fighandle, axeshandles] = myMakeMultiWaveSummaryPlot(region,data,dilateContours)
%If batchFetchWaveProp.m has been run, then the data cell array 'data' can be passed to this function for printing some summary info in the figure subplots
if nargin < 2, data = []; end
if nargin < 3 || isempty(dilateContours); dilateContours = 1; end
locationMarkers = unique(region.location);
regionsAll = splitRegion(region);

% startcounter = ones(1,numel(locationMarkers)+1);
%this counter index is used to properly slice into the optional 'data' input structure for printing relevant values on the plot.
startcounter = ones(1,numel(region.name));
if numel(locationMarkers) > 1
    for i = 2:numel(locationMarkers)
        if isempty(region.wavedata{locationMarkers(i-1)}.waveonsets)
            startcounter(locationMarkers(i)) = length(region.wavedata{locationMarkers(i-1)}.waveonsets) + 2;
        else
            startcounter(locationMarkers(i)) = length(region.wavedata{locationMarkers(i-1)}.waveonsets) + 1;
        end
    end
end


for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    try
        [fighandle, ax] = getPlotByLocation(tmpregion,locationIndex,data,startcounter,dilateContours);
        axeshandles{locationIndex}.ax = ax;
    catch exception
        if isempty(region.wavedata{locationIndex}.waveonsets)
            disp(['Region waveonsets is empty. Nothing to analyse for locationIndex = ' num2str(locationIndex)])
            startcounter(locationIndex) = startcounter(locationIndex) +1;
        else
            throw(exception);
        end
    end
end


function [fighandle, ax] = getPlotByLocation(region, locationIndex, data,startcounter,dilateContours)
fighandle = figure;
numplots = length(region.wavedata{locationIndex}.waveonsets);
cols = 6;
rows = floor(numplots/cols);
if rem(numplots,cols) > 0
    rows = rows+1;
end


[szX,szY] = size(region.image);  %assuming szY is the largest dimension
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;
end
len=length(region.contours);
f={};
for c = 1:len
    %     prg(c) = 1;
    %     figure(tfig);
    %     imagesc(prg);
    %     set(gca,'xtick',[],'ytick',[]);
    %     drawnow
    %     waitbar(c/len,wtbar);
    %     im=series1{1,1};
    ps = round(region.contours{c});
    ps=[ps(:,1) rXY*ps(:,2)];
    if rXY > 1
        idx=find(ps(:,2) == min(ps(:,2)));
        if (min(ps(:,2))-rXY) > 0
            ps(idx,2)=min(ps(:,2))-rXY;
        else
            ps(idx,2) = 1;
        end
    end
    
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    %     inp = inpolygon(subx,suby,region.contours{c}(:,1),region.contours{c}(:,2));
    %     inp = inpolygon(subx,suby,region.contours{c}(:,1),rXY*region.contours{c}(:,2));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    %     disp(num2str(c))
    f{c} = sub2ind([szX, szY],fy,fx);
end
strel_sz=[numel(unique(fx)) numel(unique(fx))];

A1=zeros([szX szY numel(region.wavedata{locationIndex}.waveonsets)]);
% i = 1; %testing
for i=1:numel(region.wavedata{locationIndex}.waveonsets)
    %     Adummy=A(:,:,region.wavedata{locationIndex}.waveonsets(i));
    %     for d=region.wavedata{locationIndex}.waveonsets(i):region.wavedata{locationIndex}.waveoffsets(i)
    % %         %         figure; imshow(A(:,:,d))
    % %         if d ~= region.wavedata{locationIndex}.waveoffsets(i)
    % %             Adummy=max(Adummy,A(:,:,d+1));
    % %         end
    %     end
    d = region.wavedata{locationIndex}.waveonsets(i):region.wavedata{locationIndex}.waveoffsets(i);
    Adummy=zeros(szX,szY);
%     figure; 
    for c=1:len
        ind = intersect(region.onsets{c},d);
%         disp(ind)
        if ~isempty(ind)
            Adummy(f{c}) = abs(ind(1)-d(end));  %set ntSpk to ntFilt if you want a contour movie that includes baseline fluctuations in between detected events.
%             imagesc(Adummy); colormap(flipud(gray));
        end
    end
%     hold off;
    %Add dummy frame to master wave list length of region.wavedata{locationIndex}.waveonsets
    A1(:,:,i)=mat2gray(Adummy);
% end  %Testing
	
	if dilateContours > 0
		%use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19px contours.
		SE = strel('rectangle', [2 2]);
		A2=imdilate(A1(:,:,i),SE);
		SE = strel('square', strel_sz(1)+3); %want close to a 3px radius to fill in the 1-2px gaps in between contours.
		A2=imclose(A2,SE);
		% figure;imshow(A2)
    else
		A2 = A1(:,:,i);
    end
    
    ax(i) = subplot(rows,cols,i);
    %   ax(i) = subaxis(rows,cols,i,'Spacing',0,'Padding',0,'Margin',0);
    imshow(A2); colormap(flipud(gray));
    title(['fr ' num2str(d(1)) ' - ' num2str(d(end))]);
    hold on
    
    %plot gray dotted outline of each region----
    for numcoords = 1:length(region.coords)
        plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
    end
    
    if i == 1
        text(0.1,0.1,region.name{locationIndex},'FontWeight','bold');
    end
    
    if isfield(region.wavedata{locationIndex},'waveprops')
        if isfield(region.wavedata{locationIndex}.waveprops,'wavedirection')
            if ~isnan(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i))
                theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) .* -1;   %origin in getWaveDirections is at lower-left corner ([y=512, x=0]). Multiplying by -1 will flip the coords to match the upper-left ([0,0]) coord system of each image suplot.
%                 theta=region.wavedata{locationIndex}.waveprops.waveorientation_radians(i);  %waveorientation from regionprops instead
%                 %theta=(region.waveprops.wavedirection.theta_radians(~isnan(region.waveprops.wavedirection.theta_radians)))+pi;
                rho=region.wavedata{locationIndex}.waveprops.wavedistance_px(i);
%                 rho=region.wavedata{locationIndex}.waveprops.wavedirection.rho_pixels(i);
                wavecentroids = round(region.wavedata{locationIndex}.wavecentroids{i}(find(~isnan(region.wavedata{locationIndex}.wavecentroids{i}(:,1))),:));
%                 rowdifference = abs(region.wavedata{locationIndex}.waveprops.wavedirection.origin_coord(1) - region.orientation.value(1));
%                 coldifference = abs(region.wavedata{locationIndex}.waveprops.wavedirection.origin_coord(2) - region.orientation.value(2));
                %         goodwavecentroids=[abs(centr(:,1)-coldifference) abs(centr(:,2)-rowdifference)];
                centr = wavecentroids(1,:);
%                 if rowdifference > centr(2)
%                 if rowdifference < centr(2)
% %                     theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) .* -1;
% %                     theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - pi/2;
% %                     theta=region.wavedata{locationIndex}.waveprops.waveorientation_radians(i) - pi/2;
%                 end
                [x,y] = pol2cart(theta,rho);
                %         centr = goodwavecentroids(1,:);
                % axes(axeshandles{locationIndex}.ax(i))
                % hold on
                plot(wavecentroids(:,1),wavecentroids(:,2),'-b')
                plot(centr(1),centr(2),'or')
                plot([centr(1) centr(1)+x],[centr(2) centr(2)+y],'-r')
                % hold off
                
                incr = szX/6;
                if ~isempty(data)
                    for ind=1:size(data,2)
                        if strcmp(data(1,ind),'waveisi.s')
                            text(0,incr,['ISI: ' num2str(cell2mat(data(startcounter(locationIndex)+i,ind))) 's'],'FontSize',9);
                        end
                        if strcmp(data(1,ind),'wavesize.um2')
                            text(0,2*incr,[num2str(round(cell2mat(data(startcounter(locationIndex)+i,ind)))) 'um2'],'FontSize',9);
                        end
                        if strcmp(data(1,ind),'wavespeed.umpersec')
                            text(0,3*incr,[num2str(round(cell2mat(data(startcounter(locationIndex)+i,ind)))) 'um/s'],'FontSize',9);
                        end
                        if strcmp(data(1,ind),'wavedir.degs')
                            text(0,4*incr,[num2str(round(cell2mat(data(startcounter(locationIndex)+i,ind)))) 'degs'],'FontSize',9);
                        end
                        if strcmp(data(1,ind),'wavedist.um')
                            text(0,5*incr,[num2str(round(cell2mat(data(startcounter(locationIndex)+i,ind)))) 'um'],'FontSize',9);
                        end
                        if strcmp(data(1,ind),'wavepathlength.um')
                            text(0,6*incr,[num2str(round(cell2mat(data(startcounter(locationIndex)+i,ind)))) 'um'],'FontSize',9);
                        end
                    end
                end
            end
        end
    end
    hold off
end
