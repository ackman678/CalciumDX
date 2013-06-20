function [Hist_all,xout,meanLatencies,nSpikes,responseArray,latencyArray]=myPETH(region,stimuliToPlot,bn,nstimuli,name,makePlots)
%myPETH - Calls getPETH after formatting and passing along a list of spike times from a 'calciumdx' or 'calciumdxcalciumdextran' formatted 'region' data structure
%region = region file.
%stimuliToPlot, stimuli numbers to plot from 'region.stimuli'
%bn = time surrounding stimuli onsets to plot for the histograms
%
%James Ackman, 2011-07-25
%
% %sdf=SDFConv(spk,AlignEv,bn,ktype,sig)
% %spk in NaN padded 2D array of spike times in msec. Each row a trial
% %AlignEv a column vector of alignment times (stimulus onset times) in msec
%
%1make AlignEv
%2cut data by idx = find between ea successive event in AlignEv
%3make spk from cut data
%4plot by stimulus, hemisphere in locationIndex
%
% TODO %  Change numel(stimuli{numStim}.stimulusParams) to a new input variable for which stimuli (like just 1st 3 before habituation) to plot

stimuli = region.stimuli;
if nargin<2 || isempty(stimuliToPlot); stimuliToPlot=1:numel(stimuli); end
if nargin<3 || isempty(bn);      bn=[-300 800];      end
if nargin<4 || isempty(nstimuli); nstimuli=1:numel(region.stimuli{stimuliToPlot(1)}.stimulusParams); end
if strcmp(nstimuli,'default'); 
nstimuliSwitch = 1;
else
nstimuliSwitch = 0;
end
if nargin<5 || isempty(name); name=''; end
if nargin<6 || isempty(makePlots); makePlots = 1; end


nt = zeros(size(region.traces));
for c = 1:size(region.traces,1)
    nt(c,:) = dfoverf(region.traces(c,:))*100;
end

if ~isempty(stimuli)
    %     numStim=1;
    for numStim = stimuliToPlot
    if nstimuliSwitch > 0; nstimuli=1:numel(region.stimuli{numStim}.stimulusParams); end
    
        if makePlots>0
        fighandle2(numStim)=figure();
        fighandle1(numStim)=figure();
        scrsize = get(0,'screensize');
        set(fighandle2(numStim),'Position',scrsize);
        set(fighandle2(numStim),'color',[1 1 1]);
        set(fighandle2(numStim),'PaperType','usletter');
        set(fighandle2(numStim),'PaperPositionMode','auto');%         numplots = numel(stimuli{numStim}.stimulusParams);
        
        set(fighandle1(numStim),'Position',scrsize);
        set(fighandle1(numStim),'color',[1 1 1]);
        set(fighandle1(numStim),'PaperType','usletter');
        set(fighandle1(numStim),'PaperPositionMode','auto');%         numplots = numel(stimuli{numStim}.stimulusParams);
        end
        numplots = numel(nstimuli);
        cols = 6;
        rows = floor(numplots/cols);
        if rem(numplots,cols) > 0
            rows = rows+1;
        end

        responseArray{numStim} = zeros(numel(region.contours),numplots);
        Hist_all = [];
%         AlignEv = [];
        meanLatencies = [];
        nSpikes = [];
%         latencyArray{numStim} = zeros(numel(region.contours),numplots);
%         for i=1:numel(stimuli{numStim}.stimulusParams)
        for i=1:numel(nstimuli)
            %         AlignEv = [AlignEv; stimuli{numStim}.stimulusParams{i}.stimulus_times(1)];
            AlignEv = repmat(stimuli{numStim}.stimulusParams{nstimuli(i)}.stimulus_times(1), size(region.traces,1),1);
            %     end
            AlignEv = AlignEv .* 1e-03;  %convert from microseconds to milliseconds
            
            spk = zeros(size(region.traces));
            for c = 1:size(spk,1)
                spk(c,region.onsets{c}) = 1;
            end
            [dummy,j] = find(spk);
            ind = find(spk);
            j = j.*region.timeres;
%            j = (j.*region.timeres) - region.timeres; %for ephys stimuli, the frames mark end of frame, not beginning.  For stimuli with <1 frame period of response latency, this can make it look like the wave form signal starts before stimulus (since each image gives the mean signal within that frame period time)
%			disp(j(1:10))  %TESTING
            spk(ind) = j;
            spk(spk == 0) = NaN;
%            figure; imagesc(spk) %TESTING
            spk = spk .* 1e03; %convert from seconds to milliseconds
            
            %     AlignEv_bkup = AlignEv;
            %     % AlignEv = AlignEv(1);
            %     AlignEv = repmat(AlignEv(1), size(spk,1),1);
            
            %[sdf,Hist_raw,xout]=getPETH(spk,AlignEv,bn,ktype,sig,binWidth,meanLatencies,nSpikes)
%             bn = [-300 800];
            [sdf,Hist_raw,xout,meanLatencies,nSpikes,cells,~,meanLatencyByCell]=getPETH(spk,AlignEv,bn,[],[],[],meanLatencies,nSpikes);
            
            if makePlots>0
            figure(fighandle1(numStim))
            ax(i) = subplot(rows,cols,i);
            bar(xout,Hist_raw,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3]);
            xlim(bn)
            axis tight
            axis square
            xlabel('Time (ms)')
            ylabel('Count')
            if i==1
                title([name 'stim ' num2str(nstimuli(i))])
            else
                title(['stim ' num2str(nstimuli(i))])
            end
            % figure; bar(xout,sdf)
            end
            
            if i == 1
                Hist_all = Hist_raw;
            else
                Hist_all = Hist_all+Hist_raw;
            end
            responseArray{numStim}(cells,i) = 1;  %for response frequency image plot in myMakeMultiPETHplot.m.  Identifies whether a cell reponded to stimulation or not. 
            
            if i== 1
                latencyArray{numStim} = meanLatencyByCell;
            else
               latencyArray{numStim} = [latencyArray{numStim} meanLatencyByCell]; 
            end
            
            %mean trace plot
            spk = ones(size(region.traces));
            [dummy,j] = find(spk);
            ind = find(spk);
            j = j.*region.timeres;
            spk(ind) = j;
            spk = spk .* 1e03; %convert from seconds to milliseconds
            hbn = bn;
            spk=spk-repmat(AlignEv,1,size(spk,2));
            [cells,frames]=find(spk>=hbn(1) & spk <= hbn(2));
            disp(i)
%            disp(frames)
            if makePlots>0
            figure(fighandle2(numStim))
            ax2(i) = subplot(rows,cols,i);
            hold on
            frames = unique(frames);
            disp(numel(frames))
            meanTrace = mean(nt(:,frames));
            plotStimuli(region,meanTrace,1,1,numStim)  %last vargin blank defaults to 0 TimeShift
            plot(frames,meanTrace,'color','k')
%             xlabs=[hbn(1):round(range(hbn)/3):hbn(2)];
%             set(gca,'xticklabel',{num2str(xlabs(1)), num2str(xlabs(2)), num2str(xlabs(4)), num2str(xlabs(4))});
            xlim([frames(1) frames(end)])
            axis square
            xlabel('Time (fr)')
            ylabel('dF/F')
            if i==1
                title([name 'stim ' num2str(nstimuli(i))])
            else
                title(['stim ' num2str(nstimuli(i))])
            end
            
            if i == 1
                meanTrace_all = meanTrace;
%                disp(meanTrace)
            else
%                disp(i)
%                disp(meanTrace)
                meanTrace_all = meanTrace_all+meanTrace;
                
            end
            end
        end
        
        
        if makePlots>0
        meanTrace_all = meanTrace_all./numel(nstimuli);
        figure(fighandle2(numStim))
        ax2(i+1) = subplot(rows,cols,i+1);
        hold on
        TimeShift = (stimuli{numStim}.stimulusParams{nstimuli(i)}.stimulus_times(1) .* 1e-03);
        frames = (frames*region.timeres*1e03) - TimeShift;
        plotStimuli(region,meanTrace_all,1,[region.timeres*1e03],numStim,TimeShift)
%         frames = ((frames-frames(1)).*region.timeres) + hbn(1);
        plot(frames,meanTrace_all,'color','k')
        xlim([frames(1) frames(end)])
%         xlabs=[hbn(1):round(range(hbn)/3):hbn(2)];
%         set(gca,'xticklabel',{num2str(xlabs(1)), num2str(xlabs(2)), num2str(xlabs(4)), num2str(xlabs(4))});
        axis square
        xlabel('Time (ms)')
        ylabel('dF/F')
        title(['mean trace all, ' region.stimuli{numStim}.description]);
        
        figure(fighandle1(numStim))
        ax(i+1) = subplot(rows,cols,i+1);
        bar(xout,Hist_all,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3]);
        title(['sum all, ' region.stimuli{numStim}.description]);
        xlim(bn)
        axis tight
        axis square
        xlabel('Time (ms)')
        ylabel('Count')
        end
    end
end

function plotStimuli(region,nt,num,timeUnits,numStim,TimeShift)
if nargin < 6 || isempty(TimeShift); TimeShift = 0; end
    if isfield(region,'stimuli')
        if ~isempty(region.stimuli)
%                 hAll = [];
%                 for numStim = 1:numel(region.stimuli)
%                     if numel(region.stimuli) > 1
%                                     mycolors = lines(numel(region.stimuli));
                        mycolors = [0.8 0.8 1.0; 0.8 1.0 0.8; 1.0 0.8 0.8; 0.6 0.6 1.0; 0.6 1.0 0.6; 1.0 0.6 0.6; 0.4 0.4 1.0; 0.4 1.0 0.4; 1.0 0.4 0.4];
%                     else
%                         mycolors = [0.8 0.8 0.8];
%                     end
                    for i=1:numel(region.stimuli{numStim}.stimulusParams)
                        x1=(region.stimuli{numStim}.stimulusParams{i}.frame_indices(1)/region.stimuli{numStim}.stimulusParams{i}.frame_times(1))*region.stimuli{numStim}.stimulusParams{i}.stimulus_times(1);
                        x2=(region.stimuli{numStim}.stimulusParams{i}.frame_indices(end)/region.stimuli{numStim}.stimulusParams{i}.frame_times(end))*region.stimuli{numStim}.stimulusParams{i}.stimulus_times(end);
                        x = [x1; x1; x2; x2]; x = (x.*timeUnits) - TimeShift;
                        y = [min(nt(num,:)); max(nt(num,:)); max(nt(num,:)); min(nt(num,:))];
                        handle_stim(i) = patch(x,y,mycolors(numStim,:));  % fill
                        %             set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
                        set(handle_stim(i),'EdgeColor',mycolors(numStim,:));  %outline
%                         set(handle_stim(i),'FaceAlpha',0.1,'EdgeAlpha',0.1)  %looks great but matlab does not export transparency well
%                         set(h1(i),'DisplayName',region.stimuli{numStim}.description{1})
                    end
%                     hAll = [hAll h1]; 
%                 end
        end
    end