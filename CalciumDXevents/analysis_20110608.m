%2011 James B. Ackman
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + pi  %fix individual wave angles
i = 12; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2)  %fix individual wave angles
i = 14; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4)  %fix individual wave angles
i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi)  %fix individual wave angles
i = 5; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2)  %fix individual wave angles
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4)  %fix individual wave angles



locationIndex=3;
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

%110323_08.mat adjustments:
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/2);
i = 7; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);
i = 9; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);
i = 10; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);
i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);
i = 7; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/8);
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);
i = 9; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi);
i = 13; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);
i = 14; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);
i = 17; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/4);
i = 20; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4);

%110308_05.mat
baseline = A_pastespecial(:,2)';
for i = 1:size(region.traces,1)
region.traces(i,:) = (region.traces(i,:) - baseline)./baseline;
end

i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/8);
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (5*pi/8);
i = 7; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);
i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4);
i = 4; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/2);
i = 5; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/4);
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 10; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/4);
i = 11; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/2);

%110613_01.mat
%manual add artifactFrames--
%20x20px ROIs, a little midline blood vessel overlap in anterior-medial couple ROIs. Some uncovered posterior-medial space in SC.L
%calciumdxdettrial.m sd=3;  sd2=1; sd3=2; nonfilt = 1; myFilter order = 2
region.artifactFrames = [257 270; 345 356; 446 465; 606 610; 638 645]

%110613_02.mat
%"imported contours" from 110613_01 to have identical coords, region.names, and ROIs
%used exact same calciumdxdettrial.m detection params as 110613_01
%manual add artifactFrames--
region.artifactFrames = [40 60; 305 324; 409 420; 446 457; 535 538; 997 1005; 1110 1130; 1360 1440]

%110308_04.mat
region.artifactFrames = [1 320; 600 1380; 1550 1670; 1742 1747; 1886 1894; 1950 2075; 2300 2600; 3050 3950; 3980 4020; 4250 4800;];
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);
i = 7; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);
i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);

%110308_06.mat
i = 4; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);
i = 7; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);
idx = find(region.location == 2); 
for i = idx  %fix SC.R, there should be no detected transients (some detected artifacts, since manual peaks wasn't used)
    region.onsets{i} = [];
    region.offsets{i} = [];
    region.transients(i) = 1;
end

%110308_07.mat
i = 2; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi);
i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);
i = 7; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);

%110304_03.mat
i = 3; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);
i = 4; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/8);

%110304_05.mat
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/2);

%110308_05.tif  small ROIs with Kalman
%setup new file with coordinates and names to those of already determined in previously analysed movie. Set current region file to region_bkup, then load file for direction extaction into workspace.
region_bkup=region;  %at the image filter selection screen, run the next few lines
clear region
%load the region.mat for importing params
region_bkup.name = region.name; region_bkup.coords = region.coords; region_bkup.wavedata = region.wavedata;
region = region_bkup; clear region_bkup;
%then run the following line at the contour filtering screen
uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 0.03],'FontSize',9,'Callback','calciumdxReadTraceParam','enable','on');

%set directions to those of already determined directions in previously analysed movie. Set current region file to region_bkup, then load file for direction extaction into workspace.
region_bkup=region;
clear region
%load the region.mat the data will be extracted from (drag and drop on matlab)
region_bkup.wavedata{2}.waveprops.wavedirection = region.wavedata{2}.waveprops.wavedirection;
region_bkup.wavedata{3}.waveprops.wavedirection = region.wavedata{3}.waveprops.wavedirection;
region = region_bkup; clear region_bkup;

%110308_06.tif   small ROIs with Kalman
region.artifactFrames(region.artifactFrames == 1453) = 1430

%110323_03.mat  small ROIs, V1.R.SC.R
i = 4; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/4);
i = 5; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (5*pi/4);
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (5*pi/4);
i = 9; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/8);
i = 10; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/2);
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 11; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);
%recorrection (some angles were still way off from visual detection, usually from wrong angle correction above) 2011-10-28
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/8);
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/8);
i = 10; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/8);
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/8);
i = 7; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/8);



%110323_05.tif   %redo with default (small) ROIs
%used same procedure as above for 110308_05.tif to setup new file with coords and names adn reimport params and directions.
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = -(1*pi/4);
%110323_08.tif   %redo with default (small) ROIs
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = -pi/4

%110323_11.tif   %small ROIs
i = 2; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/1);
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/2);
i = 10; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);
i = 11; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/4);
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/1);
i = 13; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/8);

%script for fixing previously analysed 2P files to use with batchFetchWaveprop now TSeries-12172010-1306-002_wavedet.mat
region.wavedata{1}.waveonsets = region.waveonsets; 
region.wavedata{1}.waveoffsets = region.waveoffsets;
region.wavedata{1}.wavepeaks = region.wavepeaks;
region.wavedata{1}.wavecentroids = region.wavecentroids;
region.wavedata{1}.waveprops = region.waveprops;
region.wavedata{1}.waveprops.wavedirection.origin_coord = [0 0];
region.wavedata{1}.name = region.brainarea;
region.name = {}; region.name{1} = region.brainarea;

%110613_01.tif  redo with small ROIs and medium (10px) ROIs.  10px ROIs is best
%setup new file with coordinates and names to those of already determined in previously analysed movie. Set current region file to region_bkup, then load file for direction extaction into workspace.
region_bkup=region;  %at the image filter selection screen, run the next few lines
clear region
%load the region.mat for importing params
region_bkup.name = region.name; region_bkup.coords = region.coords; region_bkup.stimuli=region.stimuli; region_bkup.artifactFrames=region.artifactFrames;
region = region_bkup; clear region_bkup;
%then run the following line at the contour filtering screen
uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 0.03],'FontSize',9,'Callback','calciumdxReadTraceParam','enable','on');

%110613_02.tif
%manually identifed artifact frames (just a few twitch movement periods) while looking at dF movie.
region.artifactFrames = [39 70; 95 120; 195 210; 300 320; 405 415; 440 465; 1110 1130; 1370 1385];

%110510_03.tif
region.artifactFrames = [199 214; 336 351; 352 357; 358 461; 735 779; 783 825; 835 875; 884 923; 932 974; 1105 1118; 1150 1159; 1425 1433];

%110510_04.tif
region.artifactFrames = [80 125; 132 175; 180 220; 420 425; 480 525; 730 775; 780 797; 798 801; 805 820; 835 850; 851 858; 859 870; 885 920; 935 970];

%must fix 110510_03.mat and 110510.mat with correct coords and contours   2011-07-21
%setup new file with coordinates and names to those of already determined in previously analysed movie. Set current region file to region_bkup, then load file for direction extaction into workspace.
region_bkup=region;  %at the image filter selection screen, run the next few lines
clear region
%load the region.mat for importing params 110510_04-2sd-0sd2-0.5sd3.mat
region_bkup.name = region.name; region_bkup.coords = region.coords;
region = region_bkup; clear region_bkup;
%then run the following line at the contour filtering screen
uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 0.03],'FontSize',9,'Callback','calciumdxReadTraceParam','enable','on');

%under calciumdxevents4 
region_bkup=region;
clear region
%load each backwards file
region_bkup.stimuli = region.stimuli; region_bkup.artifactFrames = region.artifactFrames;
region = region_bkup; clear region_bkup

%110510_02.tif synch
region.artifactFrames = [100 110; 705 725; 845 865];

%110510_01.tif unilateral
region.artifactFrames = [1 26; 64 75; 81 100; 119 126; 134 149; 161 166; 451 465; 469 481; 515 522; 523 528; 529 534; 535 540; 541 549];

%TSeries-110308.tif 2P, 2x, 20Xobj 0.134s/fr
i = 5; locationIndex = 1; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 6; locationIndex = 1; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi);  %fix individual wave angles
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (5*pi/8);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi);  %fix individual wave angles
i = 4; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (4*pi/8);  %fix individual wave angles
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);  %fix individual wave angles
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4);  %fix individual wave angles

%110323_04-fftfilt.mat
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/2);  %fix individual wave angles
i = 3; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 4; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/2);  %fix individual wave angles
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 7; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/8);  %fix individual wave angles
i = 8; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);  %fix individual wave angles
i = 9; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/4);  %fix individual wave angles


i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 5; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 6; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 9; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (7*pi/4);  %fix individual wave angles
i = 13; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);  %fix individual wave angles
i = 14; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 15; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);  %fix individual wave angles
i = 16; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (5*pi/4);  %fix individual wave angles
i = 18; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4);  %fix individual wave angles

%110309_02-fftfilt.mat
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 5; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 10; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 11; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles

%110323_01.mat
i = 3; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/4);  %fix individual wave angles
i = 4; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/8);  %fix individual wave angles

i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 3; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/2);  %fix individual wave angles
i = 6; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 10; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/8);  %fix individual wave angles
i = 11; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 13; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (5*pi/8);  %fix individual wave angles
i = 14; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/8);  %fix individual wave angles
%recorrection (some angles were still way off from visual detection, usually from wrong angle correction above) 2011-10-28
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);  %fix individual wave angles
i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/8);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/8);  %fix individual wave angles
i = 5; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/2);  %fix individual wave angles

%correction of some waveonsets for bilateral wave matching. Some of the waves start with a completely syncrhonized time. This needs to be set to the observed values (the current computer detected values has some jitter of +/- few secs.
%110323_08_defaultROIs.mat
region.wavedata{3}.waveonsets(9)=1103;
region.wavedata{2}.waveonsets(17)=2683;
region.wavedata{3}.waveonsets(19)=2683;

%110323_08_defaultROIs.mat
region.wavedata{2}.waveoffsets(14)=2349;
region.wavedata{3}.waveoffsets(16)=2349;

%correction of some waveoffsets for bilateral wave matching. Some of the waves start with a completely syncrhonized time. This needs to be set to the observed values (the current computer detected values has some jitter of +/- few secs.
%110308_05_defaultROIs_kalman2.mat
region.wavedata{2}.waveonsets(2)=760;
region.wavedata{2}.waveoffsets(2)=840;
region.wavedata{3}.waveonsets(4)=760;
region.wavedata{3}.waveoffsets(4)=840;

%/Volumes/FIRMTEK/Tim/Retinal\ Wave\ Movies/110304iT\ analysis/110304iT_5x_01.mat
region.wavedata{2}.waveonsets(2)=485;
region.wavedata{2}.waveoffsets(2)=585;
region.wavedata{3}.waveonsets(2)=485;
region.wavedata{3}.waveoffsets(2)=585;

%edit waveonset times in a couple V1-SC movies from 110323-- start times before SC, but really the waves start at same time or just after in VCtx.  But in these movies the entire SC is not neccessarily imaged at P9, whereas whole VCtx can be mostly seen depending on dye labeling extent.
%110323_04_fftfilt.mat
region.wavedata{2}.waveonsets(3) = 495;
region.wavedata{2}.waveonsets(7) = 2510;
fnm = [pathname filename];
save(fnm,'region')
%110323_01.mat
region.wavedata{2}.waveonsets(2) = 975;
region.wavedata{2}.waveonsets(4) = 1961;
fnm = [pathname filename];
save(fnm,'region')




%120221_18.mat
region.wavedata{3}.waveonsets = region.wavedata{2}.waveonsets;
region.wavedata{3}.waveoffsets = region.wavedata{2}.waveoffsets;
region.wavedata{3}.wavepeaks= region.wavedata{2}.wavepeaks;
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);  %fix individual wave angles
i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (2*pi/8);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi);  %fix individual wave angles
i = 2; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi);  %fix individual wave angles
i = 6; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 7; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 9; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 10; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/2);  %fix individual wave angles
i = 9; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (7*pi/4);  %fix individual wave angles
i = 12; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles
i = 12; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (pi/4);  %fix individual wave angles


%120221_16.mat
region.artifactFrames = [region.artifactFrames; 796 799];
i = 4; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (2*pi/4);  %fix individual wave angles
i = 3; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (2*pi/4);  %fix individual wave angles
i = 13; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/4);  %fix individual wave angles
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (5*pi/4);  %fix individual wave angles
%120221_16.mat-- 20120416, analysis revealed 3 wave angles that are way off--
i = 2; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/8);  %fix individual wave angles
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/4);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (7*pi/4);  %fix individual wave angles
i = 6; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (6*pi/4);  %fix individual wave angles


%120221_17.mat
i = 4; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (7*pi/4);  %fix individual wave angles
i = 7; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/4);  %fix individual wave angles
i = 8; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (3*pi/4);  %fix individual wave angles
i = 4; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (2*pi/4);  %fix individual wave angles
i = 6; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (1*pi/8);  %fix individual wave angles

%120221_19.mat
i = 1; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi);  %fix individual wave angles
i = 1; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (5*pi/4);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/4);  %fix individual wave angles
i = 2; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) + (3*pi/8);  %fix individual wave angles
i = 2; locationIndex = 2; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (pi/8);  %fix individual wave angles
i = 6; locationIndex = 3; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - (1*pi/8);  %fix individual wave angles