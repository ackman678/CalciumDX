Author: James B. Ackman  
Date: 2013-06-20 14:49:49  
Tags: manual, doc, protocols, methods  

# Detect cell pairs #

[fetchCorrPairs](CalciumDXevents/fetchCorrPairs.m)

The purpose of this documentation is to demonstrate how to detect significant correlations between pairs of cells/ROIs that exhibit synchronous calcium activities.

## Usage

	[region] = fetchCorrPairs(region,sig,numres,p_val,useGaussEvents,win)

default window is `win = 0.250` (250msec) if you don't use gaussEvents for gaussian smoothing.

## Example

Example using gaussian smoothing of the signal:

	[region] = fetchCorrPairs(region,1,1000,0.01,'true');

Here is an example of what will be output to the matlab command prompt:

	>> [region] = fetchCorrPairs(region,1,1000,0.01,'true');

	All cells ....................................................................................................
			Total number of cells: 397
	Percent cells in correlations: 74.3073
			Total number of pairs: 78606
		 Percent pairs correlated: 11.1315

## Explanation of outputs

The resulting returned `region` data structure will have the raw results stored at `region.userdata.corr{n}`, where *n = 1* if this is the first time you've run fetchCorrPairs, or *n = n + 1* if you've already run it *n* times.

Here are the outputs and explanations of what's returned in the `region` data structure: 

	>> region.userdata.corr{3}

	ans = 

			corr_pairs: {[8750x2 double]}	%List of cell pairs at p_val significance level
		pvalCorrMatrix: [397x397 double]	%Correlation matrix of measured p-values for all cell pairs
				params: [1x1 struct]	%parameter structure, see below

	>> region.userdata.corr{3}.params

	ans = 

		  useGaussEvents: 'true'	%Whether Gaussian smoothing of spike signal was used
			 window_sec: 1.7408		%Window in secs, approximate width of total window surrounding spike
		spkLength_frames: 5			%Windo in frames, width of each spike
				  numres: 1000		%No. of reshuffles of spike train for the Monte Carlo simulation
				   p_val: 0.0100	%Significance level set by user
					date: '20130717-083925'	%Timestamp of when the dataset was generated


To visualize the difference between *Percentage of cells in correlations* vs *Percent pairs correlated* imagine the following 6 cell network:


![][networkImg]


The cell pair list at `region.userdata.corr_pairs` along with the p-values for the pairs and their physical distances can be saved separately for doing a detailed network analysis in R or Python using the iGraph or NetworkX packages.

More on this together with a script to fetch this data structure to follow.



# Fetch Corr Pair Data #

[batchFetchCorrpairs](CalciumDXevents/batchFetchCorrpairs.m)

This batch script will output a data table for downstream plotting and analysis. 

	data = batchFetchCorrpairs(filelist)

This essentially fetches the table of percent cells correlated, number of pairs, and percent pairs correlated that is usually output to the command line from `fetchCorrPairs` as in the [Example] above.


# Plot Corr Pairs Data #

## Usage

* `myPlotCorrGraphImage(region,plotType,numLoca,alphaLevel,edgeAesthetic,datasetSelector)`
* `myPlotPvalueCorrMatrix(region, minMax, plotType, datasetSelector)`

This assumes that you've already run fetchCorrPairs.m and saved your region data file. 


## Examples


### Plot graph of corr pairs

[myPlotCorrGraphImage](CalciumDXevents/myPlotCorrGraphImage.m)

Plot a graph of nodes and edges based on the corr data:

	myPlotCorrGraphImage(region,'1',[],0.1)  %connectivity of all cells with all location types

![][img1]

	myPlotCorrGraphImage(region,'2',[],0.1)  %connectivity of all cells within one location type

![][img2]

	myPlotCorrGraphImage(region,'3',358,0.9)  %connectivity of one cell (no. 358) with all contours in all location types

![][img3]

	myPlotCorrGraphImage(region,'4',[],0.03)  %raw image with connectivity overlay

![][img4]


Save a png image of the plot (assuming you have pathname and filename defined in your workspace):

	fnm = [pathname filename];
	print('-dpng', [fnm(1:end-4) datestr(now,'yyyymmdd-HHMMSS') '.png'])


Plot edge aesthetic color coded by pvalue instead of distance:

	myPlotCorrGraphImage(region,'2',[],0.1,'pvalue',3)

![][img5]


### Plot corr matrix

[myPlotPvalueCorrMatrix](CalciumDXevents/myPlotPvalueCorrMatrix.m)

Plot a correlation matrix of the pvalues for all node pairs in your network.  With plotType = '1', the figure will have 4 plots of different min/max values for the color scale.

	myPlotPvalueCorrMatrix(region, [], '1')

![][img6]

With plotType = '1', the figure will have 4 plots of different min/max values for the color scale.  In this example the minMax values have been explicitly set to [0 0.01] instead of the default.

	myPlotPvalueCorrMatrix(region, [0 0.01], '2')
	
![][img7]


[networkImg]: assets/img/network_regular_vs_random.png "Random vs Regular networks" width="500px"
[img1]: assets/img/20130717-115851.jpg "graph distance" width="500px"
[img2]: assets/img/20130717-114721.jpg "graph distance" width="500px"
[img3]: assets/img/20130717-121446.jpg "graph distance" width="500px"
[img4]: assets/img/20130717-122440.jpg "graph distance" width="500px"
[img5]: assets/img/20130717-170425.jpg "graph pvalue" width="500px"
[img6]: assets/img/20130718-093211_figure2.png "CorrMatrix image" width="500px"
[img7]: assets/img/20130718-103109_figure2.png "CorrMatrix image" width="500px"

