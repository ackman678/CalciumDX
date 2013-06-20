%Updated by James Ackman,07/08/07
%function calciumdxprintout(fname)
% calciumdxprintout(fname)

if ~exist('region','var')
    if exist('calciumdxprefs.mat','file') == 2
        load('calciumdxprefs')
    else
        pathname = pwd;
    end
    
    [filename pathname] = uigetfile('*.mat','Select calciumdx .mat file',pathname);
    fname = fullfile(pathname,filename);
    load(fname);
    save('calciumdxprefs.mat', 'pathname','filename');
    
end
fname = [pathname filename];


spres = region.spaceres;

scrsize = get(0,'screensize');
set(gcf,'position',[scrsize(3)/2-8.5/11*0.86*scrsize(4)/2 0.07*scrsize(4) 8.5/11*0.86*scrsize(4) 0.86*scrsize(4)]);
set(gcf,'color',[1 1 1]);

f = strfind(filename,'.');
if isempty(f)
    str = filename;
else
    f = f(1);
    str = filename(1:f-1);
end

tmpres = region.timeres;

%uppper left filename------------------------------------------------------
uicontrol('Style','text','units','normalized','string',str,'position',[.1 0.95 0.43 .04],'FontSize',7,'BackgroundColor',[1 1 1]);

cl = hsv(length(region.coords));

%spectogram or periodogram of active population-------------------------------------------------
subplot('position',[0.1 0.78 0.31 0.11])
%imagesc(region.image)
hold on
% imagesize = size(region.image);
% for c = 1:length(region.contours)
%     plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',cl(region.location(c),:));
% end
% xlim([0 imagesize(2)])
% ylim([0 imagesize(1)])
% axis equal
% axis tight
% set(gca,'ydir','reverse','xtick',[],'ytick',[]);
nt = dfoverf(region.traces);
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
actvcells=find(sum(s,2)>0);
if ~isempty(actvcells)
if length(actvcells) > 1;
pwelch(mean(nt(actvcells,:)),hamming(200),[],1024,1/region.timeres)
else
pwelch(nt(actvcells,:),hamming(200),[],1024,1/region.timeres)
end    
end

% [S,F,T,P] = spectrogram(mean(nt),128,120,128,1/region.timeres);
% surf(T,F,10*log10(abs(P)),'EdgeColor','none');
% axis xy; axis tight; colormap(jet); view(0,90);
% xlabel('Time');
% ylabel('Frequency (Hz)');
% colorbar('EastOutside')
% 
% title('spectrogram(mean(region.traces))','FontWeight','bold')

% %The following line adds a scale bar
% plot([imagesize(2)-15-100/spres imagesize(2)-15],[imagesize(1)-15 imagesize(1)-15],'-k','LineWidth',2)
% xlim([0 imagesize(2)])
% ylim([0 imagesize(1)])
% axis equal
% set(gca,'ydir','reverse');
%axis off

%plot of active cells contours-----------------------------------------------
subplot('position',[0.47 0.73 0.43 0.195])
hold on
%for c = 1:length(region.coords)
%    patch(region.coords{c}(:,1),region.coords{c}(:,2),'w');
%    ct = centroid(region.coords{c});
%    tx = text(ct(1),ct(2),region.name{c});
%    set(tx,'color',cl(c,:));
%end
if isfield(region,'transients') == 1
    region.userdata.active=find(region.transients > 1)';
    region.userdata.schmutzr=find(region.transients == 2)';
    region.userdata.gdp=find(region.transients == 3)';
    region.userdata.other=find(region.transients == 4)';
    region.userdata.sink=find(region.transients == 5)';
    if ~isempty(region.userdata.active)
        for c = 1:size(region.contours,2)
            if isempty(find(region.userdata.active==c))
                plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
                
            else
                %h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'red');
                if ismember(c,region.userdata.schmutzr)
                    %set(h,'edgecolor','green');
                    h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.9 0.3 0.3]);
                elseif ismember(c,region.userdata.gdp)
                    %set(h,'edgecolor','blue');
                    h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'black');
                elseif ismember(c,region.userdata.sink)
                    %set(h,'edgecolor','blue');
                    h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'blue');
                else
                    %set(h,'edgecolor','black');
                    %         h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
                    h = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[1.0 0.5 0]);
                end
            end
        end
    else
       for c = 1:size(region.contours,2)
           plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
       end
    end
    axis equal
    imagesize = size(region.image);
    xlim([0 imagesize(2)])
    ylim([0 imagesize(1)])
    set(gca,'ydir','reverse');
    %set(gca,'color',[0 0 0]);
    set(gca,'xtick',[],'ytick',[]);
%     title('Active cells (gdp=k;sch=r)','FontSize',9,'FontWeight','bold')
    title('Active cells (other=o;sink=b)','FontSize',9,'FontWeight','bold')
    
    
    
else
%     s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
%     region.userdata.active=find(sum(s,2)>0);
    calciumdx_plot_cont(region,region.userdata.active);
end

%The following line adds a scale bar
% plot([imagesize(2)-15-100/spres imagesize(2)-15],[imagesize(1)-15 imagesize(1)-15],'-k','LineWidth',2) %100um scalebar
plot([imagesize(2)-15-50/spres imagesize(2)-15],[imagesize(1)-15 imagesize(1)-15],'-k','LineWidth',2) %50um scalebar
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
axis equal
set(gca,'ydir','reverse');
axis off

%Data table of Region, Area, #of cells, density, cell size---------------------
uicontrol('Style','text','units','normalized','string','# Of cells','position',[.01 .7 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
uicontrol('Style','text','units','normalized','string','% Active','position',[.21 .7 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
uicontrol('Style','text','units','normalized','string','Frequency (mHz)','position',[.41 .7 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
uicontrol('Style','text','units','normalized','string','Duration (sec)','position',[.61 .7 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
uicontrol('Style','text','units','normalized','string','%active synch','position',[.81 .7 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);

subplot('position',[0 .69 1 .005])
plot(xlim,[0 0],'-k')   %Draw a line
axis tight
axis off

ncells = length(region.contours);
uicontrol('Style','text','units','normalized','string',num2str(ncells),'position',[.01 .7-.03 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);


s = rast2mat(region.onsets,size(region.traces,2));

actv = sum(s,2);
pact = length(find(actv>0))/size(s,1)*100;
uicontrol('Style','text','units','normalized','string',num2str(pact,'%6.2f'),'position',[.21 .7-0.03 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
actv = actv/(size(region.traces,2)*tmpres);
actv = actv(find(actv>0)); %jba

% m_freq=num2str(mean(actv),'%6.3f'); %jba 30/07/07  in Hz
m_freq=mean(actv);
% se_freq=num2str(std(actv)/sqrt(length(actv)),'%6.3f'); %jba 30/07/07
se_freq=std(actv)/sqrt(length(actv));
actv = actv*1000; %jba make to millihertz
str = [num2str(round(mean(actv))) '±' num2str(round(std(actv)/sqrt(length(actv)))) ' (' num2str(round(min(actv))) '-' num2str(round(max(actv))) ')'];
uicontrol('Style','text','units','normalized','string',str,'position',[.41 .7-0.03 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
dur = [];
mxdur = zeros(1,length(region.contours));
mndur = zeros(1,length(region.contours));
for d = 1:length(region.contours)
    dur = [dur region.offsets{d}-region.onsets{d}];
    if ~isempty(region.offsets{d})
        mxdur(d) = max(region.offsets{d}-region.onsets{d});
        totdur(d) = sum(region.offsets{d}-region.onsets{d});
    end
end
dur = dur*tmpres;
pdur = [num2str(mean(dur),'%6.2f') '±' num2str(std(dur)/sqrt(length(dur)),'%6.2f') ' (' num2str(min(dur),'%6.1f') '-' num2str(max(dur),'%6.1f') ')'];
uicontrol('Style','text','units','normalized','string',pdur,'position',[.61 .7-0.03 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);

% m_dur=num2str(mean(dur),'%6.2f');
m_dur=mean(dur);
% se_dur=num2str(std(dur)/sqrt(length(dur)),'%6.2f');
se_dur=std(dur)/sqrt(length(dur));
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
psynch = max(sum(s,1))/size(s,1)*100;
uicontrol('Style','text','units','normalized','string',num2str(psynch,'%6.2f'),'position',[.81 .7-0.03 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);

plot(xlim,[0 0],'-k')   %Draw a line
axis tight
axis off


%----------------------------------------------------------------------------
%new info
if isfield(region,'transients') == 1
    pgdp = size(region.userdata.gdp,1)/size(region.userdata.active,1)*100;
    tpgdp = num2str(pgdp,'%6.2f');
    uicontrol('Style','text','units','normalized','string',['%gdp = ' tpgdp],'position',[.01 .7-.06 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
    
    psch = size(region.userdata.schmutzr,1)/size(region.userdata.active,1)*100;
    tpsch = num2str(psch,'%6.2f');
    uicontrol('Style','text','units','normalized','string',['%sch = ' tpsch],'position',[.21 .7-0.06 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
    
    pother = size(region.userdata.other,1)/size(region.userdata.active,1)*100;
    tpother = num2str(pother,'%6.2f');
    uicontrol('Style','text','units','normalized','string',['%other = ' tpother],'position',[.41 .7-0.06 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
    
    psink = size(region.userdata.sink,1)/size(region.userdata.active,1)*100;
    tpsink = num2str(psink,'%6.2f');
    uicontrol('Style','text','units','normalized','string',['%sink= ' tpsink],'position',[.61 .7-0.06 .18 .02],'FontSize',9,'BackgroundColor',[1 1 1]);
end

%c = c+1;
movielength = region.timeres * size(region.traces,2);

n_cells = ''; p_corr = ''; p_pairs = ''; n_pairs = ''; p_pairs_corr = '';
if isfield(region.userdata,'corr_pairs') == 1
pairs = region.userdata.corr_pairs{1};
n_cells = num2str(size(s,1));
p_corr = num2str(100*length(unique(reshape(pairs,1,numel(pairs))))/size(s,1));
n_pairs = num2str(size(s,1)*(size(s,1)-1)/2);
p_pairs_corr = num2str(100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2));
end

if isfield(region,'animaltype')
    sprintf(['filename' '\t' 'exp.type' '\t' 'genotype' '\t' 'age' '\t' 'dye' '\t' 'area' '\t' 'anesthetic' '\t' 'anesthpercent' '\t' 'field' '\t' 'zartifact' '\t' 'zdepth.um' '\t' 'frameaveraging' '\t' 'framerate.Hz' '\t' 'frameperiod' '\t' 'linesperframe' '\t' 'scanlineperiod' '\t' 'opticalzoom' '\t' 'spatialres.um' '\t' 'mov_length.sec' '\t' 'comments' '\t' 'ncells' '\t' 'pact' '\t' 'psynch' '\t' 'pgdp' '\t' 'psch' '\t' 'pother' '\t' 'psink' '\t' 'm_freq(Hz)' '\t' 'se_freq' '\t' 'm_dur(s)' '\t' 'se_dur' '\t' 'p.corr' '\t' 'n.pairs' '\t' 'pp.corr'...
        '\r' filename '\t' region.exptype '\t' region.animaltype '\t' region.age '\t' region.dye '\t' region.brainarea '\t' region.anesthetic '\t' num2str(region.anesthpercent) '\t' region.field '\t' region.zartifact '\t' num2str(region.zdepth) '\t' region.frameaveraging '\t' num2str(region.framerate) '\t' num2str(region.timeres) '\t' num2str(region.linesperframe) '\t' num2str(region.scanlineperiod) '\t' num2str(region.opticalzoom) '\t' num2str(region.spaceres) '\t' num2str(movielength) '\t' region.comments '\t' num2str(ncells) '\t' num2str(pact) '\t' num2str(psynch) '\t' num2str(pgdp) '\t' num2str(psch) '\t' num2str(pother) '\t' num2str(psink) '\t' num2str(m_freq) '\t' num2str(se_freq) '\t' num2str(m_dur) '\t' num2str(se_dur) '\t' p_corr '\t' n_pairs '\t' p_pairs_corr])
    
else
    
    
    sprintf(['filename' '\t' 'ncells' '\t' 'pact' '\t' 'psynch' '\t' 'pgdp' '\t' 'psch' '\t' 'pother' '\t' 'psink' '\t' 'm_freq(Hz)' '\t' 'se_freq' '\t' 'm_dur(s)' '\t' 'se_dur' '\t' 'p.corr' '\t' 'n.pairs' '\t' 'pp.corr'...
        '\r' filename '\t' num2str(ncells) '\t' num2str(pact) '\t' num2str(psynch) '\t' num2str(pgdp) '\t' num2str(psch) '\t' num2str(pother) '\t' num2str(psink) '\t' num2str(m_freq) '\t' num2str(se_freq) '\t' num2str(m_dur) '\t' num2str(se_dur) '\t' p_corr '\t' n_pairs '\t' p_pairs_corr])
    
end


%Network activity plots-----------------------------------------------------
%Left histogram plot of cell activity-------------------------------------
%pt = 0.65-0.03*c;
subplot('position',[0.1 0.40 0.29 0.15])
s = rast2matdur(region.onsets,region.offsets,size(region.traces,2));
x = [0 reshape(repmat(1:size(s,2)-1,2,1),1,2*size(s,2)-2) size(s,2)];
x = [x fliplr(x)];
y = reshape(repmat(sum(s)/size(s,1)*100,2,1),1,2*size(s,2));
y = [y zeros(1,2*size(s,2))];
h = patch(x,y,[0.5 0.5 0.5]);  %grey fill
set(h,'edgecolor',[0 0 0]);  %black outline
xlim([0 size(region.traces,2)])
set(gca,'xtick',[]);
title('Network activity','FontWeight','bold')
ylabel('Percent of cells active')

%Left raster plot of cell activity----------------------------------------
subplot('position',[0.1 0.25 0.29 0.15]);
hold on
for c = 1:length(region.contours)
    for d = 1:length(region.onsets{c})
        plot([region.onsets{c}(d):region.offsets{c}(d)],repmat(c,1,region.offsets{c}(d)-region.onsets{c}(d)+1),'color',[0.3 0.3 0.3]);
        %plot([region.onsets{c}(d):region.offsets{c}(d)],repmat(c,1,region.offsets{c}(d)-region.onsets{c}(d)+1),'color',cl(region.location(c),:));
    end
end
ylim([0 length(region.contours)+1]);
xlim([0 size(region.traces,2)]);
set(gca,'xtick',[]);
ylabel('Cell no.')
% set(gca,'ytick',(1:fix(length(region.contours)/100))*100)
set(gca,'ydir','reverse')
box on

%bottom left raster plot of cell fluorescence over time---------------------
% subplot('position',[0.1 0.07 0.3925 0.16]);
subplot('position',[0.1 0.07 0.40 0.16]);
%nt = dfoverf(region.traces);
% if size(nt,2) == (1000|4000);
% nt = nt(:,2:end);  %to get rid of first frame onset artifact
% imagesc(nt);
% else
imagesc(nt);
% end
%set(gca,'xtick',[]);
ylabel('Cell no.')
% set(gca,'ytick',[1 (1:fix(length(region.contours)/100))*100])
set(gca,'ydir','reverse')
set(gca,'xtick',[1 size(nt,2)]);
colorbar('EastOutside')
% mins = fix(size(nt,2)*tmpres/60);
% secs = fix(size(nt,2)*tmpres - fix(size(nt,2)*tmpres/60)*60);
secs = fix(size(nt,2)*tmpres);
% if secs == 0
%     str = '00';
% elseif secs < 10
%     str = ['0' num2str(secs)];
% else
%     str = num2str(secs);
% end
str = num2str(secs);
%box on
% set(gca,'xticklabel',{'0', [num2str(mins) ':' str]});
set(gca,'xticklabel',{'0', str});
xlabel('Time (sec)')

%Bottom example plots----------------------------------------------------------

if ~isempty(actvcells)
sdur = rast2matdur(region.onsets,region.offsets,size(nt,2));
% mxdur(find(sdur(:,size(nt,2))==1)) = 0; %old, without logical indexing
mxdur(sdur(:,size(nt,2))==1) = 0;
[mx ilong] = max(mxdur);
[mx imax] = max(sum(s,2));
indx = abs(totdur-size(nt,2)/5);
% indx(find(sum(s,2)<10))=inf; %old, without logical indexing
indx(sum(s,2)<10)=inf;
[mx iavg] = min(indx);
end

subplot('position',[0.61 0.37 0.35 0.10]);
plot(mean(nt),'color','k')
xlim([1 size(nt,2)]);
%ylim([min(nt(ilong,:))*1.05 max(nt(ilong,:))]);
set(gca,'xtick',[],'ytick',[]);
title('Mean','FontWeight','bold')
box off

if ~isempty(actvcells)
subplot('position',[0.61 0.23 0.35 0.10]);
plot(nt(ilong,:),'color','k')
xlim([1 size(nt,2)]);
ylim([min(nt(ilong,:))*1.05 max(nt(ilong,:))]);
set(gca,'xtick',[],'ytick',[]);
title('Examples','FontWeight','bold')
box off

subplot('position',[0.61 0.12 0.35 0.10]);
plot(nt(imax,:),'color','k')
xlim([1 size(nt,2)]);
ylim([min(nt(imax,:))*1.05 max(nt(imax,:))]);
set(gca,'xtick',[],'ytick',[]);
box off

subplot('position',[0.61 0.01 0.35 0.10]);
plot(nt(iavg,:),'color','k')
xlim([1 size(nt,2)]);
ylim([min(nt(iavg,:))*1.05 max(nt(iavg,:))]);
set(gca,'xtick',[],'ytick',[]);
box off
end

%uppper left print dialog-----------------------------------------------------
set(gcf,'papertype','usletter');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0.05 .05 .90 .90]);

bprint = uicontrol('Style','pushbutton','Units','normalized','String','Print!','Position',[0 0.97 .10 0.02],'FontSize',7,...
    'Callback','set(findobj(''style'',''pushbutton''),''visible'',''off''); printdlg(gcf); set(findobj(''style'',''pushbutton''),''visible'',''on'');');
fname2 = [fname(1:end-3) 'eps'];    
bsave = uicontrol('Style','pushbutton','Units','normalized','String','save','Position',[0.9 0.97 .10 0.02],'FontSize',7,...
    'Callback','set(findobj(''style'',''pushbutton''),''visible'',''off''); saveas(gcf,fname2,''epsc''); set(findobj(''style'',''pushbutton''),''visible'',''on'');');
    
% set(gcf,'closerequestfcn','clear; delete(gcf);');
set(gcf,'closerequestfcn','delete(gcf);');