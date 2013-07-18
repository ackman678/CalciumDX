function myPlotPvalueCorrMatrix(region, minMax, plotType, datasetSelector)
%myPlotCorrGraphImage(data,region,plotType,numLoca)
%Make a graph with cell outlines and physical positions as nodes and connect significant pairs with lines.
%Lines can be colorized based on cell pair distances or the actual calculated pvalues for the cell pairs.
%region must be a region data structure returned from calciumdx.m and calciumdxevents.m with corrpairs returned from fetchCorrPairs.m
%edgeAesthetic can be one of two strings: 'distance' or 'pvalue'
%EXAMPLE: myPlotCorrGraphImage([],region,'1',[],0.3, 'pvalue')  %see all connectivity with plotType 1, alpha 0.3, 'pvalue'
%EXAMPLE2: myPlotCorrGraphImage(region.userdata.corr{1}.corr_pairs{1},region,'1',[],0.3)  %see all connectivity
%EXAMPLE3: myPlotCorrGraphImage(data,region,'3',100,0.1)  %see connectivity of just cell 100
%region--
%	region data structure returned from calciumdx and calciumdxevents
%minMax-- an optional 1x2 double vector for the min and max values you want to the colormap to represent or the text string 'Auto'.  Will default to the min and max calculated pvalues for significant cell pairs.
%	[minValue maxValue]
%plotType--
%'1' = Make 4 plots with different max/min color scaling for caxis variable of colormap representing the pvalues. Auto, [minMax], [0 0.05], [0 0.01]. 
%'2' = Make single plot with the colormapping dictated by minMax
%datasetSelector--
%	

%James B. Ackman 2013-07-18 09:41:17

if nargin < 4 || isempty(datasetSelector), datasetSelector = 1; end

if nargin < 3 || isempty(plotType), plotType = '1'; end

if nargin < 2 || isempty(minMax), 
	pvalues = [];
	for i = 1:size(region.userdata.corr{datasetSelector}.corr_pairs{1},1)
			output = getPvalues(i,region.userdata.corr{datasetSelector});
			pvalues = [pvalues; output];
	end
	mn = min(pvalues); mx = max(pvalues); 
	minMax = [mn mx];
end


switch plotType    
	case '1'
		plotAll(region, minMax, datasetSelector);
	case '2'
		plotOne(region, minMax, datasetSelector);
end


function plotAll(region, minMax, datasetSelector)
figure; 
subplot(2,2,1)
imagesc(region.userdata.corr{datasetSelector}.pvalCorrMatrix)
colormap(jet)
colorbar
caxis('auto')
axis image
title('pvalue, caxis auto')

subplot(2,2,2)
imagesc(region.userdata.corr{datasetSelector}.pvalCorrMatrix)
colormap(jet)
colorbar
caxis(minMax)
axis image
title(['pvalue, caxis [' num2str(minMax) ']'])

subplot(2,2,3)
imagesc(region.userdata.corr{datasetSelector}.pvalCorrMatrix)
colormap(jet)
colorbar
caxis([0 0.05])
axis image
title('pvalue, caxis [0 0.05]')

subplot(2,2,4)
imagesc(region.userdata.corr{datasetSelector}.pvalCorrMatrix)
colormap(jet)
colorbar
caxis([0 0.01])
axis image
title('pvalue, caxis [0 0.01]')





function plotOne(region, minMax, datasetSelector);
figure; 
imagesc(region.userdata.corr{datasetSelector}.pvalCorrMatrix)
colormap(jet)
colorbar
caxis(minMax)
axis image
title(['pvalue, caxis [' num2str(minMax) ']'])






function output = getPvalues(pairNum,input)
i = input.corr_pairs{1}(pairNum,1);
j = input.corr_pairs{1}(pairNum,2);
output = input.pvalCorrMatrix(i,j);
