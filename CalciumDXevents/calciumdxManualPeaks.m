%calciumdxManualPeaks
%James Ackman, 1/6/2011
%updated 5/4/2011 to handle multiple regions
%updated 5/9/2012 the saving of region.detectorname params, to not overwrite the existing detector name and parameter information.
if ~isfield(region,'wavedata')
    region.wavedata = cell(length(region.name),1);
end
locationMarkers = unique(region.location);

for i = 1:length(region.name)
    region.wavedata{i}.name =  region.name{i};
end

spkAll = zeros(size(nt));
decAll = zeros(size(nt));
exportSpkSwitch=1;

%locationMarkers = 4;
for locationIndex = locationMarkers
    f = find(region.location == locationIndex);
    
    spk = zeros(numel(f),size(region.traces,2));
    dec = zeros(numel(f),size(region.traces,2));
    for ind = 1:numel(f)
        spk(ind,region.onsets{f(ind)}) = 1;
        dec(ind,region.offsets{f(ind)}) = 1;
    end
    
    % deltaspacing=round(30/region.timeres); %minimum 30sec Inter event interval
    if sum(sum(spk,2)) > 0
        % [hf1,trHist]=myPlotRasterHistManualPeaks(fnm,region,[],[],'true',spk,dec);
        [hf1,trHist]=myPlotRasterHistManualPeaks(fnm,region,[],[],'false',spk,dec);
        
        %{
%need either bounding box or giinput crosshairs to pick peaks
axes(hf1(2))
[xv,yv]=ginput;
pk_idx=round(xv);
        %}
        %{
onsets = cell(1,size(spk,1));
offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    onsets{c} = find(spk(c,:)==1);
    offsets{c} = find(dec(c,:)==1);
end

matAttOnOff=zeros(length(region.contours),length(region.traces));
matAttOn=zeros(length(region.contours),length(region.traces));
matAttOff=zeros(length(region.contours),length(region.traces));
numCell=length(region.contours);
for i=1:numCell
        ons=onsets{i};
        ofs=offsets{i};
    for j=1:length(ons)
%         plot([ons(j) ofs(j)-1],[i i],'Linewidth',2)
%         plot([ons(j)],[i ],'b.-')
        matAttOnOff(i,[ons(j):ofs(j)-1])=1;
        matAttOn(i,[ons(j)])=1;
        matAttOff(i,[ofs(j)])=1;
%         covMat=xcov(region.traces(3,:)',sum(matAttOnOff)/numCell'
    end
end

trHist=sum(matAttOnOff)/numCell;
        %}
        % figure; plot(trHist);
        %find peaks manually-----------------------------------
        
        %Use data brushes on figure and export values to workspace as 'waveframes' variable or 'artifactframes' variable
        % hf1 = figure;
        % hf2 = uicontrol('Position',[20 20 200 40],'String','Continue','Callback','uiresume(gcbf)');
        artifactframes=[];
        waveframes=[];
        
        %         if isfield(region.wavedata{locationIndex},'waveonsets')
        %             for wfr = 1:length(region.wavedata{locationIndex}.waveonsets)
        %                 waveframes= [waveframes region.wavedata{locationIndex}.waveonsets(wfr):region.wavedata{locationIndex}.waveoffsets(wfr)];
        %             end
        %             waveframes = waveframes';
        %         end
        
        hf2=uicontrol(hf1(1),'Style','pushbutton','Position',[20 20 200 40],'String','Continue','Callback','calciumdxManualPeaks_part2');
        hf3 = msgbox('Use Data brush while holding down shift to select multiple points that include either waveframes or artifacts, then save brushed points to workspace variable called waveframes or artifactframes','','help');
        disp('Use brushes to select data, then export as workspace variable called either waveframes or artifactframes');
        % uiwait(hf1(1));
        disp('Make sure the selected good datapoints variable is called waveframes');
        
    else
        errordlg('run auto detection first...')
        exportSpkSwitch=0;
    end
    waitfor(hf1(1))
end

if exportSpkSwitch > 0
spk = spkAll;
dec = decAll;
% [hf1,trHist]=myPlotRasterHistManualPeaks(fnm,region,[],[],'true',spk,dec);
region.detectorname=[region.detectorname ', calciumdxManualPeaks'];
end