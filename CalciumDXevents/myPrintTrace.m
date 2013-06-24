function myPrintTrace(fnm,num,nt,region,startframe,endframe,lcutoff,hipass,hipasscutoff,showStimuli)
%print trace for particular cell with subplot of high,low pass filter
%Subplot of contour data, with raw and filtered traces for the currently
%selected cell, nt(num,:), in calciumdxevents, where 'fnm' is a valid filename to save the eps figure to, 'num' is the cell number,
%'nt' is the dfoverf normalized traces, and region is the current dataset in use by
%calciumdxevents (or load 'region' independently, run dfoverf.m to generate
%'nt', and set 'num' equal to one of the cells in the dataset.
%For the filter, 'lcutoff' is between (0,1] where 1 equals the Nyquist freq. default is 0.5
%and where 'hipass' is a string of 'true or false', and 'hipasscutoff' is
%between (0,1] where 1 equals the Nyquist freq. default is 'false' and 0.005.
%example usage:
%myPrintTrace(fnm,num,nt,region,1,1000,0.2)
if nargin < 10 || isempty(showStimuli), showStimuli = 0; end
if nargin < 9 || isempty(hipasscutoff), hipasscutoff = 0.005; end
if nargin < 8 || isempty(hipass), hipass = 'false'; end
if nargin < 7 || isempty(lcutoff), lcutoff = 0.5; end
if nargin < 6 || isempty(endframe), endframe = size(nt,2); end
if nargin < 5 || isempty(startframe), startframe = 1; end

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
%Df= (0.005) * (0.5 * (1/region.timeres)) %for the hipass Desired cutoff freq

% lowWn=
% highWn=

%Adjust the Desired freq (Df) values in the following two lines for desired
%filter. b = fir1(n,Wn,'ftype'), where: n = order of window (10 is the
%default order of the low pass filter, but may need to be incr to 20 or
%higher for certain applications), Wn is a number between 0 and 1, where 1 corresponds to the Nyquist frequency 
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-Nyq may
%be needed for low pass filter.
% lcutoff = 0.5;
lcutofffreq = lcutoff*Nyq;
% hipass = 'false';
if strcmp(hipass,'true')
xf=filtfilt(fir1(300,hipasscutoff,'high'),1,xf);  %need just for baseline correction typically
end
xf=filtfilt(fir1(10,lcutoff,'low'),1,xf);  %lcutoff close to 1 (Nyquist frequency) will be less lowpass filtering, practically the same as no low pass filter. Lower number (below 0.1) will be more robust filtering

subplot(3,1,2);
% plot((1:length(nt))*tres,nt(num,:),'color','k');  %why did this one work with length before?
hold on
if showStimuli > 0
   plotStimuli(region,nt,num,1) 
end
plot((startframe:endframe),nt(num,startframe:endframe),'color','k');
xlim([startframe endframe]) 
xlabel('Frame no.'); 
ylabel('dF/F'); 
title(['cell ' num2str(num) ', raw trace']);

subplot(3,1,3);
% plot((1:length(nt))*tres,xf,'color','k'); %why did this one work with length before?
% plot((startframe:endframe)*tres,xf(startframe:endframe),'color','k'); 
hold on
if showStimuli > 0
   plotStimuli(region,xf,1,tres) 
end
plot((startframe:endframe)*tres,xf(startframe:endframe),'color','k');
xlim([startframe endframe]*tres)
% xlim([startframe endframe])
% set(gca,'XTickLabel',([startframe:endframe]')*tres)
xlabel('Time (s)');
ylabel('dF/F'); 
title(['low pass filter cutoff = ' num2str(lcutoff) '*Nyquist freq (' num2str(lcutofffreq) '), hipass filter =' hipass]);

figname = ['cell' num2str(num)]; fname2 = [fnm(1:end-4) figname 'fr' num2str(startframe) '-' num2str(endframe)]; saveas(gcf,fname2,'epsc');

function plotStimuli(region,nt,num,timeUnits)
    if isfield(region,'stimuli')
        if ~isempty(region.stimuli)
%                 hAll = [];
                for numStim = 1:numel(region.stimuli)
%                     if numel(region.stimuli) > 1
%                                     mycolors = lines(numel(region.stimuli));
                        mycolors = [0.8 0.8 1.0; 0.8 1.0 0.8; 1.0 0.8 0.8; 0.6 0.6 1.0; 0.6 1.0 0.6; 1.0 0.6 0.6; 0.4 0.4 1.0; 0.4 1.0 0.4; 1.0 0.4 0.4];
%                     else
%                         mycolors = [0.8 0.8 0.8];
%                     end
                    for i=1:numel(region.stimuli{numStim}.stimulusParams)
                        x1=(region.stimuli{numStim}.stimulusParams{i}.frame_indices(1)/region.stimuli{numStim}.stimulusParams{i}.frame_times(1))*region.stimuli{numStim}.stimulusParams{i}.stimulus_times(1);
                        x2=(region.stimuli{numStim}.stimulusParams{i}.frame_indices(end)/region.stimuli{numStim}.stimulusParams{i}.frame_times(end))*region.stimuli{numStim}.stimulusParams{i}.stimulus_times(end);
                        x = [x1; x1; x2; x2]; x = x.*timeUnits;
                        y = [min(nt(num,:)); max(nt(num,:)); max(nt(num,:)); min(nt(num,:))];
                        handle_stim(i) = patch(x,y,mycolors(numStim,:));  % fill
                        %             set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
                        set(handle_stim(i),'EdgeColor','none');  %outline
%                         set(handle_stim(i),'FaceAlpha',0.1,'EdgeAlpha',0.1)  %looks great but matlab does not export transparency well
%                         set(h1(i),'DisplayName',region.stimuli{numStim}.description{1})
                    end
%                     hAll = [hAll h1]; 
                end
        end
    end