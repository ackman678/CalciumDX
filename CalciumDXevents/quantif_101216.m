%quantif_101216
%Testing /100913/TSeries-09132010-1150-040
[filename, pathname] = uigetfile({'*.tif'}, 'Choose image to open');
% if ~isstr(filename{k})
if ~isstr(filename)
    return
end
% fnm = [pathname filename{k}];
fnm = [pathname filename];

data=bfopen(fnm);
numPlanes = size(data{1, 1},1);
sz = size(data{1, 1}{1,1});


tic
d1=zeros([sz numPlanes],'uint16');
% d1=zeros([sz numPlanes]);
for i=1:numPlanes
   d1(:,:,i)=data{1, 1}{i,1};
end
toc
clear data

%without fuzzy frame delete (just auto artifact detect, and auto event detect (calciumdxeventdettrial.m) for cell 102--  4/8 poss wave transients are detected. Problem is with onset calc from preevent baseline determination

%compare diff of rawdata with diff of 2nd order Hann filtered data (the filter order i'm using for some of the event detection algorithm)
cellnum=102;
tr = myfilter(region.traces(cellnum,:),2);
figure; subplot(4,1,1)
plot(region.traces(cellnum,:))
title('rawtrace')
subplot(4,1,2)
plot(tr)
title('myfilter, hanning, order 2')
subplot(4,1,3)
plot(diff(region.traces(cellnum,:)))
title('diff rawtrace')
subplot(4,1,4)
plot(diff(tr))
title('diff filtered trace')
%----------------------------

%----now test 2D median filt-----
figure; imagesc(d1(:,:,887))

L = medfilt2(d1(:,:,887));
figure; imagesc(L);
title('3x3 default medfilt2 int16')

tic
d2=zeros([sz numPlanes],'uint16');
for i=1:numPlanes
   d2(:,:,i)=medfilt2(d1(:,:,i));
end
toc
implay(mat2gray(d2))
%---------------------------------


%----Now duplicate the region .mat file for /100913/Tseries_040 and perform trace reading on median filtered array, then rerun calciumdxevents
%dupl necessary data and settings from old region file, and clear the trace and event data
region_dupl=region;
clear region;
region=region_dupl;
clear region_dupl;
region.traces=[];
region.onsets=[];
region.offsets=[];
region.transients=[];
region.artifactFrames=[];
region.onsets = cell(1,length(region.contours));
region.offsets = cell(1,length(region.contours));

%read in file
%Testing /100913/TSeries-09132010-1150-040
[filename, pathname] = uigetfile({'*.tif'}, 'Choose image to open');
% if ~isstr(filename{k})
if ~isstr(filename)
    return
end
% fnm = [pathname filename{k}];
fnm = [pathname filename];

data=bfopen(fnm);
numPlanes = size(data{1, 1},1);
sz = size(data{1, 1}{1,1});

series1=zeros([sz numPlanes],'uint16');
% series1=zeros([sz numPlanes]);
for i=1:numPlanes
   series1(:,:,i)=medfilt2(data{1, 1}{i,1});  %use default medfilt2 from above (which has a 3x3 default). Median filter is the best one for getting rid of detector shot noise.
end
clear data

oldfolder=cd('./TraceReaders')
% [tr, trhalo, param] = calciumdxTR_ReadTracesPrTIFFmedianfilt([],region,series1);
[tr, trhalo, param] = calciumdxTR_ReadTracesPrTIFF([],region,series1);
region.traces = tr;
region.transients = ones(1,size(region.traces,1));
cd(oldfolder)
%This worked pretty well,, improved S/N for event detection considerably.
%--------------------------------------------------------------------------

%---Now read in good averaged image, and compare with rest, wave, and movement artifact frames using xcorr
A=imread('/Users/ackman/Figures and images/2photon/100913/TSeries-09132010-1150-040/AVG_TSeries-09132010-1150-040_F0(fr874-881).tif');
figure; imagesc(A); title('good avg image')

A1=d1(:,:,887); %select a wave frame from medfilt2 raw data
figure; imagesc(A1); title('fr887, def medfilt2')

A2 = medfilt2(d1(:,:,887));  %refilter the def medfilt2 wave frame. A 3x3 or 6x6 Gaussian will be better
figure; imagesc(A2);
title('medfilt2(fr887 def medfilt2)')

A3 = medfilt2(A); %medfilt2 on the good avg image from above to blur and smooth. A 3x3 or 6x6 Gaussian will be better
figure; imagesc(A3);
title('good avg image with 3x3 default medfilt2 int16')

m=mean(reshape(A,1,numel(A)));  %Get the mean value and std from the Good Avg image
sd=std(reshape(A,1,numel(A)));
rn = m + sd.*randn(256,512);  %Make random noise image with same mean and std and same size as our data
fig; imagesc(rn); title('randn noise with mean and std of avg image A')

B=xcorr2(A,rn);
figure; imagesc(B); title('A vs randn')


B=normxcorr2(double(A1),double(A1));  %normalized xcorr2. corr values between [-1 1]. This is an autocorr. Many values above 0.5, and close to or at 1. Highest values directly in image center (256,512). Rest of values offset to upperleft quadrant, different than xcorr2
figure; imagesc(B); title('A1 vs A1, normxcorr2')
figure; hist(reshape(B,1,numel(B))); title('A1 vs A1, normxcorr2')
[y1,x1]=find(B==max(B(:)))

A4=padarray(A1, [512 512]);  %gives same visual result, with this much smaller template
B=normxcorr2(A1,A4);
figure; imagesc(B); title('A1 padded vs A1')

B=normxcorr2(A1,rn);  %normalized xcorr2 of wave frame vs rand normal noise.  No values above 0.4, and highest values not at center.
figure; imagesc(B); title('A1 vs rn, normxcorr2')
figure; hist(reshape(B,1,numel(B))); title('A1 vs rn, normxcorr2')
[y1,x1]=find(B==max(B(:)))

B=normxcorr2(A1,A);  %normalized xcorr2 of wave frame vs rand normal noise.  No values above 0.4, and highest values not at center.
figure; imagesc(B); title('A1 vs A, normxcorr2')
figure; hist(reshape(B,1,numel(B))); title('A1 vs A, normxcorr2')
[y1,x1]=find(B==max(B(:)))

A2=d1(:,:,57); %a wave frame from earlier in the movie, before some zdrift happenened, since the avg image A is from just before fr887.
figure; imagesc(A2); title('fr57, def medfilt2')
B=normxcorr2(A2,A);  %normalized xcorr2. Maybe the normxcorr vals aren't terribly informative. as the max value of A1 vs rn is higher than A2 vs the avg image. But highest values do cluster near center of image. So the best similarity measure is the center of mass in the image. 
figure; imagesc(B); title('A2, fr57waveframe vs A, normxcorr2')
figure; hist(reshape(B,1,numel(B))); title('A2, fr57 waveframe vs A, normxcorr2')
[y1,x1]=find(B==max(B(:)))

%--back to xcorr2
B=xcorr2(A1,rn);  %normalized xcorr2 of wave frame vs rand normal noise.  No values above 0.4, and highest values not at center.
figure; imagesc(B); title('A1 vs rn, xcorr2')
figure; hist(reshape(B,1,numel(B))); title('A1 vs rn, xcorr2')
[y1,x1]=find(B==max(B(:)))

B=xcorr2(A1,A);  %normalized xcorr2 of wave frame vs rand normal noise.  No values above 0.4, and highest values not at center.
figure; imagesc(B); title('A1 vs A, xcorr2')
figure; hist(reshape(B,1,numel(B))); title('A1 vs A, xcorr2')
[y1,x1]=find(B==max(B(:)))

A2=d1(:,:,57); %a wave frame from earlier in the movie, before some zdrift happenened, since the avg image A is from just before fr887.
B=xcorr2(A2,A);  %normalized xcorr2. Maybe the normxcorr vals aren't terribly informative. as the max value of A1 vs rn is higher than A2 vs the avg image. But highest values do cluster near center of image. So the best similarity measure is the center of mass in the image. 
figure; imagesc(B); title('A2, fr57waveframe vs A, xcorr2')
figure; hist(reshape(B,1,numel(B))); title('A2, fr57 waveframe vs A, xcorr2')
[y1,x1]=find(B==max(B(:)))


%compare padded array with non-padded array results-- No difference. Histo looks different due to all the zero values, but the main corr image array and max value and location is unchanged.
A4=padarray(A1, [512 512]);
B=normxcorr2(A,A4);
figure; imagesc(B); title('A1 padded vs A')
figure; hist(reshape(B,1,numel(B))); title('A1 padded, fr887 waveframe vs A, normxcorr2')
[y1,x1]=find(B==max(B(:)))
max(B(:))

B=normxcorr2(A,A1);
figure; imagesc(B); title('A1 vs A')
figure; hist(reshape(B,1,numel(B))); title('A1, fr887 waveframe vs A, normxcorr2')
[y1,x1]=find(B==max(B(:)))
max(B(:))

% now try with xcorr2 comparison-- Same here. No difference. Histo looks different due to all the zero values, but the main corr image array and max value and location is unchanged.
A4=padarray(A1, [512 512]);
B=xcorr2(A,A4);
figure; imagesc(B); title('A1 padded vs A')
figure; hist(reshape(B,1,numel(B))); title('A1 padded, fr887 waveframe vs A, xcorr2')
[y1,x1]=find(B==max(B(:)))
max(B(:))

B=xcorr2(A,A1);
figure; imagesc(B); title('A1 vs A')
figure; hist(reshape(B,1,numel(B))); title('A1, fr887 waveframe vs A, xcorr2')
[y1,x1]=find(B==max(B(:)))
max(B(:))

%---Gaussian blur code-----------------------
%{
w=fspecial('gaussian',[r c],sig); %Gaussian lowpass filter of size rxc and standard deviation sig. The defualts are 3x3 and 0.5. A single number instead of [r c] specifies a square filter.
IF=imfilter(A,w,'replicate');
figure; imshow(IF)
%}

%---the following is equal to a Gaussian blur with sigma = 6 in ImageJ
w=fspecial('gaussian',[51 51],6);
IF=imfilter(A,w,'replicate');
figure; imagesc(IF); title('w=fspecial("gaussian",[51 51],6));')
figure; imagesc(A); title('raw data')


accuracy=2e-4
kernelradius=ceil(3*sqrt(-2*log(accuracy)))+1
ImageJGaussBlurKernelSidelength=2*kernelradius-1; %for Gaussian Blur with sigma=6, this kernel is 51x51 (sidelength=51)

%use more than one avg image (up to 3) and loop through all frames using ea as a template. Keep highest corr values and keep frames that match one of the 3 templates
%----------------------------------------------
%------image xcorr on gaussian blurred data----
%{
%this is the d1 from opened from above with median noise reduction
d1=zeros([sz numPlanes],'uint16');
% d1=zeros([sz numPlanes]);
for i=1:numPlanes
   d1(:,:,i)=medfilt2(data{1, 1}{i,1});  %use default medfilt2 from above (which has a 3x3 default). Median filter is the best one for getting rid of detector shot noise.
end
%}

% d2=zeros([sz numPlanes],'uint16');
sz=size(d1);
d2=zeros(sz);
% w=fspecial('gaussian',[51 51],6);
w=fspecial('gaussian',[27 27],3);
% d1=zeros([sz numPlanes]);
for i=1:numPlanes
   d2(:,:,i)=imfilter(d1(:,:,i),w);  %use default medfilt2 from above (which has a 3x3 default). Median filter is the best one for getting rid of detector shot noise.
end

Afilt=imfilter(A,w);

tic
corrvals=zeros(1000,3);
for i = 1:numPlanes
    %     B=normxcorr2(A,d2(:,:,i));
    B=normxcorr2(Afilt,d2(:,:,i));
    %     B=xcorr2(Afilt,d2(:,:,i));
    [corrvals(i,1),corrvals(i,2)]=find(B==max(B(:)));
    corrvals(i,3)=max(B(:));
end
corrvals1=corrvals;
toc

%try procrustes distance measure analysis
tic
corrvals=zeros(1000,1);
for i = 1:numPlanes
    d=procrustes(Afilt,d2(:,:,i));
    corrvals(i,1)=d;
end
corrvals2=corrvals;
toc

figure;
subplot(2,1,1)
plot(corrvals1(:,3)); title('normxcorr2 max corr value'); grid on
subplot(2,1,2)
plot(corrvals2(:,1)); title('procrustes dissimilarity measure'); grid on

figure; plot(corrvals(:,3)) %corr values are near 0.8 for matched frames, especially during non-wave frames. But the big incr in intensities during wave frames causes it to be less matched, so the corr val actually decr during wave frames.
figure; plot(corrvals(:,1)) %this one is the x translation in max corr value (anterior-posterior direction-- max displacement during movements). This one might match might match the wave frames pretty well actually.

tic; B=normxcorr2(Afilt,d2(:,:,887)); toc; %1.090808sec
tic; B=xcorr2(Afilt,d2(:,:,887)); toc;   %5.872145sec
d=procrustes(Afilt,d2(:,:,887));


%--normxcorr2 loop for 3 different templates-------------------------------
A1=imread('/Users/ackman/Figures and images/2photon/100913/TSeries-09132010-1150-040/AVG_TSeries-09132010-1150-040_F0(fr8-10)_medfilt1gauss3.tif');
A2=imread('/Users/ackman/Figures and images/2photon/100913/TSeries-09132010-1150-040/AVG_TSeries-09132010-1150-040_F0(fr257-264)_medfilt1gauss3.tif');
A3=imread('/Users/ackman/Figures and images/2photon/100913/TSeries-09132010-1150-040/AVG_TSeries-09132010-1150-040_F0(fr874-881)_medfilt1gauss3.tif');

Afilt=A1;
corrvals=zeros(1000,3);
for i = 1:size(d2,3)
    B=normxcorr2(Afilt,d2(:,:,i));
    [corrvals(i,1),corrvals(i,2)]=find(B==max(B(:)));
    corrvals(i,3)=max(B(:));
end
corrvals1=corrvals;

Afilt=A2;
corrvals=zeros(1000,3);
for i = 1:size(d2,3)
    B=normxcorr2(Afilt,d2(:,:,i));
    [corrvals(i,1),corrvals(i,2)]=find(B==max(B(:)));
    corrvals(i,3)=max(B(:));
end
corrvals2=corrvals;

Afilt=A3;
corrvals=zeros(1000,3);
for i = 1:size(d2,3)
    B=normxcorr2(Afilt,d2(:,:,i));
    [corrvals(i,1),corrvals(i,2)]=find(B==max(B(:)));
    corrvals(i,3)=max(B(:));
end
corrvals3=corrvals;

%--this simple fourieralign function only takes 28sec, but doesn't do a very good job---------
tic
Afilt=A3;
corrvals4=zeros(1000,2);
for i=1:size(d2,3)
shift=fourieralign(Afilt,d2(:,:,i));
corrvals4(i,1)=shift(1);
corrvals4(i,2)=shift(2);
end
toc


%x-values
figure;
subplot(3,1,1)
plot(corrvals1(:,1)); title('normxcorr2 template fr8-10'); grid on
subplot(3,1,2)
plot(corrvals2(:,1)); title('normxcorr2 template fr257-264'); grid on
subplot(3,1,3)
plot(corrvals3(:,1)); title('normxcorr2 template fr874-881'); grid on

corrvals=corrvals3;
myPlotRasterHist(fnm,region,[],[],'true',corrvals)


%y-values
figure;
subplot(3,1,1)
plot(corrvals1(:,2)); title('normxcorr2 template fr8-10'); grid on
subplot(3,1,2)
plot(corrvals2(:,2)); title('normxcorr2 template fr257-264'); grid on
subplot(3,1,3)
plot(corrvals3(:,2)); title('normxcorr2 template fr874-881'); grid on

%corr values
figure;
subplot(3,1,1)
plot(corrvals1(:,3)); title('normxcorr2 template fr8-10'); grid on
subplot(3,1,2)
plot(corrvals2(:,3)); title('normxcorr2 template fr257-264'); grid on
subplot(3,1,3)
plot(corrvals3(:,3)); title('normxcorr2 template fr874-881'); grid on

%distance(hypotenuse)
figure;
subplot(3,1,1)
plot(sqrt(corrvals1(:,1).^2 + corrvals1(:,2).^2)); title('normxcorr2 template fr8-10'); grid on
subplot(3,1,2)
plot(sqrt(corrvals2(:,1).^2 + corrvals2(:,2).^2)); title('normxcorr2 template fr257-264'); grid on
subplot(3,1,3)
plot(sqrt(corrvals3(:,1).^2 + corrvals3(:,2).^2)); title('normxcorr2 template fr874-881'); grid on

%----merge waveforms for template1 and 2
corrvals4=corrvals2(:,1)+corrvals3(:,1);
myPlotRasterHist(fnm,region,[],[],'true',corrvals4)

%---use waveform input to delete event onsets for artifacts
h=myPlotRasterHistFrameCorr(fnm,region,[],[],'true',corrvals3)
axes(h(1))
[x,y] = ginput(1)
badframes=find(corrvals3(:,1) < y);

for i=1:length(region.contours)
    [c,idx] = setdiff(region.onsets{i}, badframes);
    region.onsets{i}=c;
    region.offsets{i}=region.offsets{i}(idx);
    if isempty(region.onsets{i})
        region.transients(i)=1;
    end
end

g=myPlotRasterHistFrameCorr(fnm,region,[],[],'true',corrvals3)

m1 = mean(region.traces,1);
% m2 = repmat(m1,size(region.traces,1),1);
% imagesc(dfoverf(region.traces - m2));
figure; plot(m1) %use linear fit to determine slope
y=-0.31x+610;
slope=0.31;

y=(1:1000).*-0.31 + 610;


% m2 = repmat(m1.*slope,size(region.traces,1),1);
m2=repmat(y,size(region.traces,1),1);

m3=region.traces - m2;
figure; plot(mean(m3,1))
nt = dfoverf(region.traces);
figure; imagesc(region.traces)
figure; imagesc(m3)
nt = dfoverf(m3);
figure; imagesc(nt);
set(gca,'clim',[-100 100])


x=m1;
% x=region.traces;
ntFilt2=[];
windowSize=100;
for i=1:size(x,1)
     y=filter(ones(1,windowSize)/windowSize,1,x(i,:)'); %same as rolling window mean filter of length 'windowSize'
    ntFilt2=[ntFilt2; y'];
end

ntFilt2=smooth(x,50);

figure; 
% subplot(2,1,1)
plot(m1)
hold on
plot(ntFilt2,'color','r')
% subplot(2,1,2)
% plot(ntFilt2)

cftool(m1)  %save fitted exponential model to workspace
coeffvalues(fittedmodel1)


%print out some traces for fig1--------------------------------------------
%100913/TSeries-09132010-1150-040_imagearraymedfilt2.mat
%repeat for cell 29,443,55,77, 83, 148,152,166,220,239
myPrintTrace(fnm,29,nt,region,850,910)


%--build script to find true Fraction of ROIs involved in ea wave
numel(find(spk(:,850:900) > 0)) / length(region.contours)     %this gives a value of 0.4951. Compare with peak of 0.4756 in histogram
%0.5016 vs ~0.40, 0.16 vs 0.10, 0.09 vs 0.09, 0.11 vs 0.09, 0.34 vs 0.20, 0.14 vs 0.11, 0.02 vs 0.02, 0.21 vs 0.19
myPlotRasterHist(fnm,region,[],[],[])

%--test out baseline filtering---------------------------------------------
%myPrintTrace(fnm,num,nt,region,1,1000,0.2)  
m1 = mean(region.traces,1);
myPrintTrace(fnm,1,m1,region,1,length(region.traces),0.5,'true')

%--Build script to intersect nt trace image values with events-------------
%----must first do hipass filtering for baseline correction of region.traces to improve the intersected nt trace image of events (otherwise the colormap is scaled to the first waves, if there are strong baseline changes in mean fluorescence)---
Nyq=0.5*(1/region.timeres);
hipasscutoff = 0.005;
ntFilt=zeros(size(region.traces));

if 900 > size(region.traces,2)
    highfilterorder=round((size(region.traces,2)/3)) - 1;  %data must of length of greater than 3x filter order.
elseif size(region.traces,2) == 900
    highfilterorder=round((size(region.traces,2)/3)) - 2;
else
    highfilterorder = 300;
end

for i=1:size(region.traces,1)
   xf=filtfilt(fir1(highfilterorder,hipasscutoff,'high'),1,region.traces(i,:));
   ntFilt(i,:)=xf;
end
ntFilt=mat2gray(ntFilt);


%------Now intersect the new ntFilt values with event times for a thresholded event matrix coded by event amplitudes
figure; imshow(mat2gray(ntFilt)); colormap(jet)
ntSpk=zeros(size(nt));
for c=1:size(nt,1)
    if ~isempty(region.onsets{c})
        for d=1:length(region.onsets{c})
            %            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = region.traces(c,region.onsets{c}(d):region.offsets{c}(d));
            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = ntFilt(c,region.onsets{c}(d):region.offsets{c}(d));
        end
    end
end
figure; imshow(mat2gray(ntSpk)); colormap(flipud(gray))   %this flips the colormap so we have a white background
figure; imshow(mat2gray(ntSpk)); colormap(jet)

%2011-02-13
%print out some traces for fig2--------------------------------------------
%/101008/TSeries-10082010-1458-003_rectmeshgrid_wavedet2.mat
%repeat for cell 52 78 104 124 142 162 263
%fr120-170
myPrintTrace(fnm,104,nt,region,120,170)
