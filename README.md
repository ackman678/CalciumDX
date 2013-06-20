Date: 2013-06-19 15:53:40  
Author: James B. Ackman


#Image processing, ROI detection, and calcium trace reading

These are the instructions on how to use `calciumdx.m` to perform the initial reading in of your image time-series, make ROIs, and read raw fluorescence trace data.

1. At matlab command prompt change into the calciumdx folder then type `calciumdx`
1. Click 'Open'. For CCD recording select your .tif file. For 2P, select first .tif image in your series. Wait for image series to open, may take up to a minute to complete.
2. Select to average frames together for base image. Want at least 5-10 consecutive frames that are stable (no xyz movements) and are largely non-active (not much activity). May need to open series as virtual stack in imageJ to know which frame range to select
	* Select the anterior medial reference point in the image.
3. Observe the spatial resolution and temporal resolution values to make sure they are correct. For CCD recording you must manually enter the values (i.e. 2.27275 um/px for 5x obj, or 4.5455 um/px for 2.5x obj and 0.2 sec/frame)
4. You can now adjust brightness and contrast so you can see edges of your labeled region in your image.
5. Click 'Add' button to add each region to the image using polygon outline. (i.e outline each labeled hemisphere or brain region independently). Then 'left click' to add points to surrount the structure, when complete 'right click' to complete outline. Repeat 'Add' button for each structure (i.e. 1 for ea hemisphere).
6. Click 'Next'
7. Enter names for each region ('SC.R', 'SC.L', 'V1.R', etc). Region 1 will be everything outside your labeled regions, this could be simply named 'craniotomy'.
8. Click 'Next'
9. Select 'calciumdxIF_none_' from drop down menu and click 'Filter'
10. Click and reclick the 'Hide' button.
11. Click the bottom'Next' to continue on (unless performing automatic cell detection and manual ROI drawing inside labeled regions).
12. When asked to draw rectangular grid, select 'Yes'.
13. When asked to draw rectanglar ROI grid for the 1st region (probably the 'craniotomy' background region), select 'No'.
14. When asked to draw rectangular ROI grid for each labele d region, select 'Yes'.
15. Enter a ROI size of 20x20px for CCD recordings. Leave at the default value for 2Pphoton based recordings. Close the resulting figure image windows.
16. On the main application window select 'ReadTracesPrTIFF' and click 'Next>>'
17. Select '_None_' from the Signal detector drop down menu and click 'Detect!'.
18. Click 'Finish' and save the resulting .mat file.

#Calcium event detection

For this we will use `calciumdxevents4.m`

1. Click 'Artifact rem w/FFT' button.
2. Click 'Detect Artifacts' button.
3. Click 'Detect all' button.
4. Click 'Save' to save file.
5. Click 'Manual Peaks'. Follow instructions using the brush tool and export variable either called 'waveframes' or 'artifactframes'. This more important for 2P recorded movies where artifacts may be more evident, or movies with low frequency of real activity where peaks of artifacts could be picked up as waves in downstream analysis. If the primary majority of detected transients are real wave activity then just click 'Continue' without brushing and exporting any variables. This will just accept all frames as possible frames containing real activity. Then the next couple scripts will be used for sorting out the real wave activity containing periods.
6. Click 'Save' to save file again.


#Calcium wave detection

1. Copy each of the following lines in sequence and paste to the matlab command prompt:

region = calciumdxDetectWaves(region);
region = calciumdxDetectWavesRefine(region);
region = getWaveCentroids(region);
region = getWaveSizeDistance(region);
region = getActiveFraction(region);
region = getWaveSpeeds(region);
region = getWaveDirections(region);
data = batchFetchWaveProp({filename},region); 
[fighandle, axeshandles] = myMakeMultiWaveSummaryPlot(region,data);
	
%or batch run lines 3-9 with following one liner:

	region = getWaveCentroids(region); region = getWaveSizeDistance(region); region = getActiveFraction(region); region = getWaveSpeeds(region); region = getWaveDirections(region); data = batchFetchWaveProp({filename},region); [fighandle, axeshandles] = myMakeMultiWaveSummaryPlot(region,data);

2. Now fix any necessary wave angles by looking at the WaveSummaryPlot and entering the following command:

	%fix individual wave angles.  i = wave no. to fix; locationIndex = region location (usually 2 for SC.R or 3 for SC.L); enter the subtraction or addition angle at end of line in radians in the range `[+2pi -2pi] (i.e. +/-(pi), +/-(pi/2), +/-(pi/4), +/-(3*pi/4), +/-(3*pi/2)`.
	i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi)  
	data = batchFetchWaveProp({filename},region); myMakeMultiWaveSummaryPlot(region,data); %copy output to excel

3. Save .mat file and the resulting multiwaveplots:
	save(fnm,'region')
	fname2 = [fnm(1:end-4) 'multiwaveplotRhemi' '.eps']; printfig('epsc',fname2)
	fname2 = [fnm(1:end-4) 'multiwaveplotLhemi' '.eps']; printfig('epsc',fname2)



# Summary Plots

	myPlotRasterHist([],region,[],[],[],'true','true'); fname2 = [fnm(1:end-4) 'rastHistSmooth' '.eps']; printfig('epsc',fname2)
	myPlotRasterHist([],region,[],[],[],'false','true'); fname2 = [fnm(1:end-4) 'rastHist' '.eps']; printfig('epsc',fname2)
	calciumdxprintout4.m

Wave dir plot------------------------------------------
	%print out and save wave direction plots with following command sequence (copy and paste to command line):
	locationIndex=3;  %change this no. to the hemisphere you want to plot
	theta=region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(~isnan(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians));
	rho=region.wavedata{locationIndex}.waveprops.wavedistance_px(~isnan(region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians));
	rho = rho .* region.spaceres;
	
	%should add the following lines (need to be edited) to flip coord system (for display purposes) of the plots so that 0,90,180,270 are pointing towards post,lateral,anterior,medial directions for both hemsipheres respectively. Use same lines as from batchFetchWaveProp.m and myMakeMultiWaveSummaryPlot.m
	%         centr = goodwavecentroids(1,:);
	%         disp(['centr = ' num2str(centr)])
	%         disp(['rowdiff = ' num2str(rowdifference)])
	%         if results.wavedirection.origin_coord(1) > centr(2)
	%               theta = theta .* -1;
	%         end
	
	[x,y] = pol2cart(theta,rho);
	figure; compass(x,y)
	title([region.name(locationIndex) ' all waves, wavedistance in um'])
	fname2 = [fnm(1:end-4) 'wavedirLhemi' '.eps']; printfig('epsc',fname2)
	%-----------------END wave dir plot---------------------------------------