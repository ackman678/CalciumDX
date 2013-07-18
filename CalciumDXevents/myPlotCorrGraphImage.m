function myPlotCorrGraphImage(region,plotType,numLoca,alphaLevel,edgeAesthetic,datasetSelector)
%myPlotCorrGraphImage(data,region,plotType,numLoca)
%Make a graph with cell outlines and physical positions as nodes and connect significant pairs with lines.
%Lines can be colorized based on cell pair distances or the actual calculated pvalues for the cell pairs.
%region must be a region data structure returned from calciumdx.m and calciumdxevents.m with corrpairs returned from fetchCorrPairs.m
%EXAMPLE: myPlotCorrGraphImage(region,'1',[],0.3, 'pvalue')  %see all connectivity with plotType 1, alpha 0.3, 'pvalue'
%EXAMPLE2: myPlotCorrGraphImage(region,'1',[],0.3)  %see all connectivity
%EXAMPLE3: myPlotCorrGraphImage(region,'3',100,0.1)  %see connectivity of just cell 100
%region--
%	region data structure returned from calciumdx and calciumdxevents.  The relevant data should be stored within this structure at 'region.userdata.corr{n}.corr_pairs', where n is 1 or the datasetSelector value.
%plotType--
%'1' = plot all connectivity between all contours in all locations,
%'2' = plot connectivity between all contours within one location, uses numLoca to determine the region.location
%'3' = plot connectivity between one contour with all contours in all locations, uses numLoca to determine the contour
%'4' = plot region.image with all connections
%numLoca--
%	region.location number to plot for plot type '2' or a cell number (that gives a unique region.location number) for plot type '3'.  Defaults to 2, (assumes your first region.location 1 is nothing, like craniotomy, neuropil, or areas in slice with no labeled cells)
%alphaLevel--
%	alpha transparency level for plotting the edges. Defaults to 0.5.
%edgeAesthetic--
%	%edgeAesthetic can be one of two strings: 'distance' or 'pvalue'. What you want the edge coloring to represent. Defaults to 'distance'
%datasetSelector--
%	A single number n, representing what region.userdata.corr{n} dataset you want to plot. Defaults to 1.
%James B. Ackman 2013-01-05 14:46


%=====Setup default parameters=========

if nargin < 7 || isempty(datasetSelector), datasetSelector = 1; end

if nargin < 6 || isempty(edgeAesthetic), edgeAesthetic = 'distance'; end

if nargin < 5 || isempty(alphaLevel), alphaLevel = 0.5; end

if nargin < 4 || isempty(numLoca), numLoca = 2; end

if nargin < 3 || isempty(plotType), plotType = '1'; end

%if nargin < 1 || isempty(data), data = region.userdata.corr{datasetSelector}.corr_pairs{1}; end
data = region.userdata.corr{datasetSelector}.corr_pairs{1};

%The following is important for getting the distances right if the data pixel dimensions are not equivalent
%And below the scripts will assume wherever 'rXY' is used, that it is szX (m dimension) which must be scaled up.
%the following assumes that the modulus of raster scanned data is 0 (equally divisible image size) and that for CCD images the ratio of image dimensions is either equivalent or not equally divisible
[szX,szY] = size(region.image);  %assuming szY is the largest dimension and szX may or may not need to be scaled.
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;
end

%===Contour Graph plot of connectivity edges over image of cell contours as nodes=========

switch plotType    
	case '1'
		[dist,edgeList,pvalues] = plotAll(data,region,numLoca,rXY,datasetSelector);
	case '2'
		[dist,edgeList,pvalues] = plotOneLocation(data,region,numLoca,rXY,datasetSelector);
	case '3'
		[dist,edgeList,pvalues] = plotOneCell(data,region,numLoca,rXY,datasetSelector);
	case '4'
		[dist,edgeList,pvalues] = plotRawImage(data,region,numLoca,rXY,datasetSelector);
end

%===Plot the edges and colorize based on edgeAesthetic==========

mycolors = jet(numel(dist));  %make colormap based on number of pairs to plot returned from the plot commands above
if strcmp(edgeAesthetic, 'distance')
	edgeData = [edgeList dist];  %make Nx3 list of pairs and their edgeAesthetic (physical distance or p-value)
	mx = max(dist); mn = min(dist);
elseif strcmp(edgeAesthetic, 'pvalue')
	edgeData = [edgeList pvalues];  %make Nx3 list of pairs and their edgeAesthetic (physical distance or p-value)
	mx = max(pvalues); mn = min(pvalues);
end

edgeData = sortrows(edgeData,3);   %sort the Nx3 list of pairs based on their edgeAesthetic
h2=[];
for i = 1:size(edgeData,1)  %Plot the pairs in order based on their sorted edgeAesthetic and color their connections with the colormap
    centr1 = centroid(region.contours{edgeData(i,1)});  %Get centroid location for cell 1
    centr2 = centroid(region.contours{edgeData(i,2)});  %Get centroid location for cell 2
    h2(i)=patch([centr1(1) centr2(1)],[centr1(2)*rXY centr2(2)*rXY],[0 0 0]);  %Draw patch line object connecting the cells
    set(h2(i),'EdgeColor',mycolors(i,:))
    %     set(h2(i),'EdgeColor','k')
    set(h2(i),'FaceAlpha',alphaLevel,'EdgeAlpha',alphaLevel)
end
colorbar
caxis([mn mx])
title(['plotType' num2str(plotType) ',alpha' num2str(alphaLevel) ';' edgeAesthetic ',[min,max]=[' num2str(mn) ',' num2str(mx) ']'])
%printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])






%==========Plotting functions=============================================================

function [dist,edgeList,pvalues] = plotAll(data,region,numLoca,rXY,datasetSelector)

%the following sets which region.location you want to analyse
% numLoca = unique(region.location);
%         numLoca = 2;
figure;
set(gcf,'color',[1 1 1]);
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
h1=[];
for c = 1:size(region.contours,2)
	%             if region.location(c) == numLoca
	%                 plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
	%                 %    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
	%                 %    set(h1(c),'FaceAlpha',0.1,'EdgeAlpha',0.1)
	%             else
	plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
	%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
	%    set(h1(c),'FaceAlpha',0,'EdgeAlpha',0)
	%             end
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)]);
ylim([0 imagesize(1)*rXY]);
set(gca,'ydir','reverse');
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
axis off

% i=4;

dist=[];
edgeList=[];
pvalues = [];
% i = 1;
distAll =[];
for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	distAll = [distAll; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
end

for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	%             if region.location(data(i,1)) == numLoca && region.location(data(i,2)) == numLoca
	edgeList = [edgeList; data(i,:)];	
	if isfield(region.userdata.corr{datasetSelector},'pvalCorrMatrix')
		output = getPvalues(i,region.userdata.corr{datasetSelector});
		pvalues = [pvalues; output];
	end
	dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
	%             end
end







function [dist,edgeList,pvalues] = plotOneLocation(data,region,numLoca,rXY,datasetSelector)

%the following sets which region.location you want to analyse
% numLoca = unique(region.location);
%         numLoca = 2;
figure;
set(gcf,'color',[1 1 1]);
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
h1=[];
for c = 1:size(region.contours,2)
	if region.location(c) == numLoca
		plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
		%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
		%    set(h1(c),'FaceAlpha',0.1,'EdgeAlpha',0.1)
	else
		plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[1 1 1]);
		%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
		%    set(h1(c),'FaceAlpha',0,'EdgeAlpha',0)
	end
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)]);
ylim([0 imagesize(1)*rXY]);
set(gca,'ydir','reverse');
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
axis off

% i=4;

dist=[];
edgeList=[];
pvalues = [];
% i = 1;
distAll =[];
for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	distAll = [distAll; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
end

for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	if region.location(data(i,1)) == numLoca && region.location(data(i,2)) == numLoca
		edgeList = [edgeList; data(i,:)];
		if isfield(region.userdata.corr{datasetSelector},'pvalCorrMatrix')
			output = getPvalues(i,region.userdata.corr{datasetSelector});
			pvalues = [pvalues; output];
		end
		dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
	end
end







function [dist,edgeList,pvalues] = plotOneCell(data,region,numLoca,rXY,datasetSelector)
%the following sets which region.location you want to analyse
% numLoca = unique(region.location);
%         numLoca = 2;
figure;

set(gcf,'color',[1 1 1]);
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
h1=[];
for c = 1:size(region.contours,2)
	if c == numLoca
		plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0 0 0]);
		%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
		%    set(h1(c),'FaceAlpha',0.1,'EdgeAlpha',0.1)
	else
		plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
		%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
		%    set(h1(c),'FaceAlpha',0,'EdgeAlpha',0)
	end
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)]);
ylim([0 imagesize(1)*rXY]);
set(gca,'ydir','reverse');
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
axis off

% i=4;

dist=[];
edgeList=[];
pvalues=[];
% i = 1;
distAll =[];
for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	distAll = [distAll; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
end

for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	if data(i,1) == numLoca || data(i,2) == numLoca
		edgeList = [edgeList; data(i,:)];
		if isfield(region.userdata.corr{datasetSelector},'pvalCorrMatrix')
			output = getPvalues(i,region.userdata.corr{datasetSelector});
			pvalues = [pvalues; output];
		end
		dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
	end
end




function [dist,edgeList,pvalues] = plotRawImage(data,region,numLoca,rXY,datasetSelector)

%the following sets which region.location you want to analyse
% numLoca = unique(region.location);
%         numLoca = 2;
figure;
imshow(region.image,[])
%{
set(gcf,'color',[1 1 1]);
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
h1=[];
for c = 1:size(region.contours,2)
%             if region.location(c) == numLoca
%                 plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
%                 %    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
%                 %    set(h1(c),'FaceAlpha',0.1,'EdgeAlpha',0.1)
%             else
		plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
		%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
		%    set(h1(c),'FaceAlpha',0,'EdgeAlpha',0)
%             end
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)]);
ylim([0 imagesize(1)*rXY]);
set(gca,'ydir','reverse');
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
axis off
%}


% i=4;

dist=[];
edgeList=[];
pvalues = [];
% i = 1;
distAll =[];
for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	distAll = [distAll; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
end

for i = 1:size(data,1)
	centr1 = centroid(region.contours{data(i,1)});
	centr2 = centroid(region.contours{data(i,2)});
	%             if region.location(data(i,1)) == numLoca && region.location(data(i,2)) == numLoca
	edgeList = [edgeList; data(i,:)];
	if isfield(region.userdata.corr{datasetSelector},'pvalCorrMatrix')
		output = getPvalues(i,region.userdata.corr{datasetSelector});
		pvalues = [pvalues; output];
	end
	dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
	%             end
end




function output = getPvalues(pairNum,input)
i = input.corr_pairs{1}(pairNum,1);
j = input.corr_pairs{1}(pairNum,2);
output = input.pvalCorrMatrix(i,j);
