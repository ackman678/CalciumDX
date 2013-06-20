%the following script will shift the onsets values to the left by one frame, useful when calciumdxdettrial is used for auto detection and sets the onsets too far forward in time.

for i = 1:length(tr);
if spk(num,i) == 1;
spk(num,i) = 0;
spk(num,i-1) = 1;
end
end

%the following will remove the first frame (recording onset artifact), for visualization purposes
nt = nt(:,2:end);

plot(filter(hann(5),1,nt(num,:)),'r') %matlab hann filter
plot(myfilter(nt(num,:),2),'g')   %rosa's hann filter

%make calciumdxevents2 window active and overplot with filtered trace
plot(myfilter(nt(num,:),2),'black')   %rosa's hann filter


%LOAD .mat file
[file path] = uigetfile();
fname = [path file];
load(fname);

%SAVE file
save(fname,'region','-v6')

%Could do [file path] = uigetfile(); fname(1) = [path file]; to make fname file list for selection of a bunch of .mat files for further input into a batch analysis program



%example for importing data from another .mat file
[tfilename, tpathname] = uigetfile(); 
tfnm=fullfile(tpathname,tfilename); 
tmp = load(tfnm); 
region.contours=tmp.region.contours; 
region.name=tmp.region.name; 
region.coords=tmp.region.coords; 

%Example code from hevSave for utilizing, the onset/offset array variables
region.onsets = cell(1,size(spk,1));
region.offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    region.onsets{c} = find(spk(c,:)==1);
    region.offsets{c} = find(dec(c,:)==1);
end

%For fixing the length of region.onsets, offsets  if first frame has been removed.
for c = 1:size(region.onsets,2)
region.onsets{c} = region.onsets{c} - 1;
region.offsets{c} = region.offsets{c} -1;
end

%get frequency data for only synchronized/oscillating networks over ages (GDP selected movies)

%network synchrony/activation pattern by spatial distribution distance 

%FIND GDP cells (visual inspection of network synch peaks)
region.userdata.gdp=find_gdp(region);
s=rast2matdur(region.onsets,region.offsets,size(region.traces,2));
region.userdata.active=find(sum(s,2)>0);
region.userdata.drugs=('nbqxapv');
region.userdata.age=('10');
save ValJu031A1116cont region

%HISTOGRAM-onsets set up onsets in a binary matrix for histograms/rasterplots.  This one sets up the binary matrix taking in account just the onset timepoints.
sz=1000;
s = zeros(length(region.onsets),sz);
for c = 1:size(s,1)
    for f = region.onsets{c}
        s(c,f) = 1;
    end
end

nice_hist(s);

%HISTOGRAM-duration that includes the transient duration times(ea. timepoint = 1).
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
x = [0 reshape(repmat(1:size(s,2)-1,2,1),1,2*size(s,2)-2) size(s,2)];
x = [x fliplr(x)];
y = reshape(repmat(sum(s)/size(s,1)*100,2,1),1,2*size(s,2));
y = [y zeros(1,2*size(s,2))];
h = patch(x,y,[0 0 0.5]);


%compute number of active cells in movie
x = sum(s,2)
x(x > 0)
numel(x(x > 0))

%compute number of synchronized active cells at one time
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
max(sum(s,1))/size(s,1)
nice_hist(s);
xlim([0 size(region.traces,2)])

%save graphics object
fname2 = [fname(1:end-3) 'eps']; 
saveas(gcf,fname2,'epsc'); %must be epsc to save color

%print a formatted report, for easy copy and paste into excel cells
sprintf([ncells '\t' pact '\t' psynch '\t' pgdp '\t' psch '\t' pother])

%SELECT cells in contour plot to identify cell indices
selectcell(region,[])

%print trace for a particular cell
%num = 47;  %optional if calciumdxevents and cell of interest is already selected
%watch out cause right now its written to delete first frame. 
tres = region.timeres; figure(); plot((2:length(nt))*tres,nt(num,2:end),'color','k'); xlabel('Time (s)'); ylabel('dF/F'); Title(['cell ' num2str(num)]);
figname = ['cell' num2str(num)]; fname2 = [fnm(1:end-4) figname '.' 'eps']; saveas(gcf,fname2,'epsc');





%print trace for particular cell with subplot of high,low pass filter--------------------------------------
%Subplot of contour data-----------------------------------------------------
figure();
subplot(3,1,1);
hold on;
for c = 1:size(region.contours,2)
        plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
end
h = patch(region.contours{num}([1:end 1],1),region.contours{num}([1:end 1],2),'black');
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
%set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);


%figure();
xf=nt(num,:);
tres = region.timeres;
Nyq=0.5*(1/region.timeres);
%Wn = Df/(0.5 * Sf);  %Df = desired cutoff freq; Sf = sampling frequency (1/region.timeres)

%Adjust the Desired freq (Df) values in the following two lines for desired
%filter. b = fir1(n,Wn,'ftype'), where: n = order of window (10 is the
%default order of the low pass filter, but may need to be incr to 20 or
%higher for certain applications), Wn is a number between 0 and 1, where 1 corresponds to the Nyquist frequency 
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-Nyq may
%be needed for low pass filter.
lcutoff = 0.5;
lcutofffreq = lcutoff*Nyq;
hipass = 'false';
if strcmp(hipass,'true')
xf=filtfilt(fir1(300,0.005/Nyq,'high'),1,xf);  %need just for baseline correction typically
end
xf=filtfilt(fir1(10,(lcutoff*Nyq)/Nyq,'low'),1,xf);  %1 or 1.5 works good

subplot(3,1,2);
plot((1:length(nt))*tres,nt(num,:),'color','k');
xlabel('Time (s)'); 
ylabel('dF/F'); 
Title(['cell ' num2str(num) ', raw trace']);

subplot(3,1,3);
plot((1:length(nt))*tres,xf,'color','k');
%xlabel('Time (s)'); 
ylabel('dF/F'); 
Title(['low pass filter cutoff = ' num2str(lcutoff) '*Nyquist (' num2str(lcutofffreq) '), hipass =' hipass]);

figname = ['cell' num2str(num)]; fname2 = [fnm(1:end-4) figname '.' 'eps']; saveas(gcf,fname2,'epsc');
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------




%Print drug comparison traces---------------------------------------------------------------------------
%print trace for particular cell with subplot of high,low pass filter--------------------------------------
%Subplot of contour data-----------------------------------------------------
figure();
%-----1
num=num1;
nt=nt1;
region=region1;
subplot(3,2,1);
hold on;
for c = 1:size(region.contours,2)
        plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
end
h = patch(region.contours{num}([1:end 1],1),region.contours{num}([1:end 1],2),'black');
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
%set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);

%------3
%figure();
xf=nt(num,:);
tres = region.timeres;
Nyq=0.5*(1/region.timeres);
%Wn = Df/(0.5 * Sf);  %Df = desired cutoff freq; Sf = sampling frequency (1/region.timeres)

%Adjust the Desired freq (Df) values in the following two lines for desired filter
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-2 may be needed for low pass filter. 10 is the default order of the low pass filter, but may need to be incr to 20 or higher for certain applications.
%xf=filtfilt(fir1(300,0.005/Nyq,'high'),1,xf);  %need just for baseline correction typically
xf=filtfilt(fir1(10,1/Nyq,'low'),1,xf);  %1 or 1.5 works good

subplot(3,2,3);
plot((1:length(nt))*tres,xf,'color','k');
%xlabel('Time (s)'); 
ylabel('dF/F'); 
title(['cell ' num2str(num)]);

%------5
subplot(3,2,5);
plot((1:length(nt))*tres,nt(num,:),'color','k');
xlabel('Time (s)'); 
ylabel('dF/F'); 

%-----2
num=num1;
nt=nt2;
region=region2;
subplot(3,2,2);
hold on;
for c = 1:size(region.contours,2)
        plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
end
h = patch(region.contours{num}([1:end 1],1),region.contours{num}([1:end 1],2),'black');
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
%set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);

%------4
%figure();
xf=nt(num,:);
tres = region.timeres;
Nyq=0.5*(1/region.timeres);
%Wn = Df/(0.5 * Sf);  %Df = desired cutoff freq; Sf = sampling frequency (1/region.timeres)

%Adjust the Desired freq (Df) values in the following two lines for desired filter
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-2 may be needed for low pass filter. 10 is the default order of the low pass filter, but may need to be incr to 20 or higher for certain applications.
%xf=filtfilt(fir1(300,0.005/Nyq,'high'),1,xf);  %need just for baseline correction typically
xf=filtfilt(fir1(10,1/Nyq,'low'),1,xf);  %1 or 1.5 works good

subplot(3,2,4);
plot((1:length(nt))*tres,xf,'color','k');
%xlabel('Time (s)'); 
ylabel('dF/F'); 
title(['cell ' num2str(num)]);

%------6
subplot(3,2,6);
plot((1:length(nt))*tres,nt(num,:),'color','k');
xlabel('Time (s)'); 
ylabel('dF/F'); 


figname = ['cell' num2str(num)]; fname2 = [fnm(1:end-4) figname '.' 'eps']; saveas(gcf,fname2,'epsc');
%-------------------------------------------------------------------------