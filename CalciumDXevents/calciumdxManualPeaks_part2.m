%calciumdxManualPeaks_part2
%James Ackman, 1/20/2011

delete(hf2);
if ishandle(hf3)
    close(hf3)
end

%======START optional lines================================================
%the following optional several lines if you want to pass waveframe or artifactframe indices to additonal regions in your movie (like if activity or artifacts are happening at the same times in multiple regions)

% if locationIndex == 2
% %    waveframes = waveframes2;
%     waveframes = waveframes1;
% end
%if locationIndex == 3
%%%    waveframes = [waveframes2; waveframes];
%%%    artifactframes = artifactframes2;
%    artifactframes = artifactframes1;
%end
%if locationIndex == 4
%%%    artifactframes = artifactframes3;
%	artifactframes = artifactframes1;
%end
%=======END optional lines=================================================


if isempty(waveframes) && isempty(artifactframes)
    button = questdlg('There are no waves... Delete all detected transients?','No waveframes selected, waveframes variable empty');
    if strcmp(button,'Yes')
        d=1:length(region.traces);
%         close(hf1(1))
%     elseif strcmp(button,'No')
    else
        d=[];
%         close(hf1(1))
%     else
%         hf2=uicontrol(hf1(1),'Style','pushbutton','Position',[20 20 200 40],'String','Continue','Callback','uiresume(gcbf)');
%         uiwait(hf1(1));
%         delete(hf2);
%         if isempty(waveframes)
%             button = questdlg('Delete all detected transients in region.onsets?','Still no waveframes selected, waveframes variable empty');
%             if strcmp(button,'Yes')
%                 d=1:length(region.traces);
%                 close(hf1(1))
%             else
%                 d=[];
%                 close(hf1(1))
%             end
%         else
%             close(hf1(1))
%         end
    end
end

if ~isempty(artifactframes) && isempty(waveframes)    
    waveframes = setxor(artifactframes(:,1), [1:length(region.traces)])';
    c=intersect(find(trHist > 0),waveframes(:,1)); %this will limit the selected possible waveframes, to just those that could have occured during a wave (those frames with 1 or more active ROIs).
    d = setxor(c, [1:length(region.traces)]);
if ishandle(hf1(1))
%     close(hf1(1))
end
    trHist(d)=0;
    figure;
    a(1) = subplot(2,1,1);
    plot(trHist)
    a(2) = subplot(2,1,2);
    imagesc(nt)
    linkaxes(a,'x')
    zoom xon
end


if isempty(artifactframes) && ~isempty(waveframes)
    c=intersect(find(trHist > 0),waveframes(:,1)); %this will limit the selected possible waveframes, to just those that could have occured during a wave (those frames with 1 or more active ROIs).
    d = setxor(c, [1:length(region.traces)]);
if ishandle(hf1(1))
%     close(hf1(1))
end
    trHist(d)=0;
    figure;
    a(1) = subplot(2,1,1);
    plot(trHist)
    a(2) = subplot(2,1,2);
    imagesc(nt)
    linkaxes(a,'x')
    zoom xon
end

%----need to integrate f here---
onsets = cell(1,size(spk,1));
offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    onsets{c} = find(spk(c,:)==1);
    offsets{c} = find(dec(c,:)==1);
end

for i=1:size(spk,1)
    [c,idx] = setdiff(onsets{i}, d);
    onsets{i}=c;
    offsets{i}=offsets{i}(idx);
    if isempty(onsets{i})
        region.transients(f(i))=1;
    end
end


%concatenate with f here
for c = 1:size(spk,1)
    spkAll(f(c),onsets{c}) = 1;
    decAll(f(c),offsets{c}) = 1;
end

% for i=1:length(region.contours)
%     if isempty(onsets{i})
%         region.transients(i)=1;
%     end
% end

% [hf1,trHist]=myPlotRasterHistManualPeaks(fnm,region,[],[],'true',spk,dec);
% region.detectorname='unsupervised calciumdxdettrial and calciumdxManualPeaks';
region.wavedata{locationIndex}.waveframes = waveframes;
region.wavedata{locationIndex}.artifactframes = artifactframes;

%{
%similar to manualpeak.m
axes(h(1))
[x,y] = ginput(1)
%}

%{
%from hevZoom2.m
point1 = get(gca,'CurrentPoint');
finalRect = rbbox;
point2 = get(gca,'CurrentPoint');

candx = sort([point1(1) point2(1)]);
if range(candx) > 1
    xlimits = sort([point1(1) point2(1)]);
else
    xlimits = get(gca,'UserData');
end
xlim(xlimits)
%}
close(hf1(1))
