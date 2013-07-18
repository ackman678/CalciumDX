function myPlotPvalueCorrMatrix(region, minMax, plotType, datasetSelector)
%myPlotPvalueCorrMatrix Plot pvalue correlation matrix.
%myPlotPvalueCorrMatrix(region, minMax, plotType, datasetSelector)
%
%Examples:
%myPlotPvalueCorrMatrix(region, [], '1')  %use default minMax which will be the range of pvalues for **significant** cell pairs, use 1st plot type which provides 4 variations. Default dataset will be used at region.userdata.corr{1}
%myPlotPvalueCorrMatrix(region, [0 0.01], '2', 3)  %manually set minMax, use 2nd plot type, and used 3rd dataset at region.userdata.corr{3}
%
%**USE**
%%region must be a region data structure returned from calciumdx.m and calciumdxevents.m with corrpairs returned from fetchCorrPairs.m
%Options:
%region-- region data structure returned from calciumdx and calciumdxevents and fetchCorrPairs
%minMax-- an optional 1x2 double vector for the min and max values [minValue maxValue] you want to the colormap to represent or the text string 'Auto'.  Will default to the min and max calculated pvalues for significant cell pairs.
%plotType--
%'1' = Make 4 plots with different max/min color scaling for caxis variable of colormap representing the pvalues. Auto, [minMax], [0 0.05], [0 0.01]. 
%'2' = Make single plot with the colormapping dictated by minMax
%datasetSelector-- A single number n, representing what region.userdata.corr{n} dataset you want to plot. Defaults to 1.
%	
%Versions:
%James B. Ackman 2013-07-18 09:41:17
%See also:
%fetchCorrPairs, myPlotPvalueCorrMatrix, batchFetchCorrPairs

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
