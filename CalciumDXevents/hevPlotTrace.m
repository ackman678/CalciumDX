% delete(h1); delete(h2); delete(h3); delete(h4); delete(h5);
% clear trmenu cmenu coffmenu h1 h2 h3 h4 h5
% set(ax(2),'xlim',xlimits)
% set(ax(3),'xlim',xlimits)
% % axes(ax(3))
% 
% hold off

button_state = get(btSameYAxis,'Value');
if button_state == get(btSameYAxis,'Max')
	% Toggle button is pressed, take appropriate action
	ylimits = [minY maxY];
	set(ax(3),'YLim',ylimits);
elseif button_state == get(btSameYAxis,'Min')
	% Toggle button is not pressed, take appropriate action
	ylimits = [min(nt(num,:)) max(nt(num,:))];
	set(ax(3),'YLim',ylimits);
end


[xlimits, ~]=queryXYLimits(ax(3));

set(cnt,'facecolor',[0 0 0]);
set(cnt(num),'facecolor',1-cl(region.location(num),:));

idx = region.transients(1,num);
set(rtransients,'value',idx);

% set(radio(1),'value',0);
% set(radio(2),'value',0);
% set(radio(3),'value',0);
% set(radio(4),'value',0);
% set(radio(idx),'value',1);

%set(ax(1),'xlim',xlimits)
%set(ax(2),'xlim',xlimits)
%set(fig,'CurrentAxes',ax(3))

hold off;
trmenu = uicontextmenu;
h1 = uimenu(trmenu, 'Label', 'Add event', 'Callback', 'hevAddEvent');

% %plot(nt(num,:),'uicontextmenu',trmenu) %the default plot
% %-----START rmBaseline
% m1 = mean(region.traces,1);
% m2 = repmat(m1,size(region.traces,1),1);
% nt = dfoverf(region.traces - m2);
% %-----END rmBaseline
% 


%-------Next two lines are for overplotting with filtered trace------
%-------added by JBA Aug 1, 2008---------------------------
% h8 = plot(filtfilt(fir1(5,0.2,'low'),1,nt(num,:)),'k'); %overplot lowpass filtered trace in black
% hold on
% h9 = line([1 len],[5 5],'LineStyle','--');
% h10 = line([1 len],[-5 -5],'LineStyle','--');
plot(ax(3),nt(num,:),'color',[.75 .75 .75],'uicontextmenu',trmenu); %default with gray trace


% -------Plot the onsets and offsets of events-------------------------
hold on
f = find(spk(num,:)==1);
g = find(dec(num,:)==1);
cmenu = zeros(1,length(f));
coffmenu = zeros(1,length(f));
for c = 1:length(f)
    cmenu(c) = uicontextmenu;
    h6 = uimenu(cmenu(c), 'Label', 'Delete event', 'Callback', ['selev = ' num2str(c) '; hevDeleteEvent;']);
    h7 = uimenu(cmenu(c), 'Label', 'Move onset', 'Callback', ['selev = ' num2str(c) '; hevMoveEvent;']);
    h11 = uimenu(cmenu(c), 'Label', 'Fuzzy delete all for this frame', 'Callback', ['selev = ' num2str(c) '; hevFuzzyDeleteAll;']);
    if c > 1
        uimenu(cmenu(c), 'Label', 'Combine with previous', 'Callback', ['selev = ' num2str(c) '; hevCombineEvent;']);
    end
    h3 = plot(f(c),nt(num,f(c)),'or','uicontextmenu',cmenu(c));
    
    coffmenu(c) = uicontextmenu;
    h4 = uimenu(coffmenu(c), 'Label', 'Move offset', 'Callback', ['selev = ' num2str(c) '; hevMoveEventOff;']);
    h5 = plot(g(c),nt(num,g(c)),'og','uicontextmenu',coffmenu(c));
end

if showArtifacts==1
if isfield(region,'artifactFrames')
    if ~isempty(region.artifactFrames)
     plot(region.artifactFrames(:,1),nt(region.artifactFrames(:,1)),'+k');
    plot(region.artifactFrames(:,2),nt(region.artifactFrames(:,2)),'*k');
    end
end
end

if showStimuli==1
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
                        x = [x1; x1; x2; x2];
                        y = [min(nt(num,:)); max(nt(num,:)); max(nt(num,:)); min(nt(num,:))];
                        handle_stim(i) = patch(x,y,mycolors(numStim,:));  % fill
                        %             set(h1(i),'EdgeColor',mycolors(numStim,:));  %outline
%                         set(handle_stim(i),'EdgeColor','none');  %outline
                        set(handle_stim(i),'EdgeColor',mycolors(numStim,:));  %outline
                        set(handle_stim(i),'FaceAlpha',0.3,'EdgeAlpha',0.3)  %looks great but matlab does not export transparency well
%                         set(h1(i),'DisplayName',region.stimuli{numStim}.description{1})
                    end
%                     hAll = [hAll h1]; 
                end
        end
    end
end

set(ax(3),'XLim',xlimits);
set(ax(3),'YLim',ylimits);


%---------------------------------------------------------------------
% set(gcf,'KeyPressFcn','hevButtonDown')
% set(gca,'buttondownfcn','hevZoom')
% set(fig,'KeyPressFcn','hevButtonDown','Interruptible','off','BusyAction','cancel')
