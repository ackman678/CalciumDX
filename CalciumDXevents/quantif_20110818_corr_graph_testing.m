%quantif_20110818_corr_graph_testing.m
%get correlations and plot contour graph of connectivity and corr dist histograms
[corr_pairs,useGaussEvents,win,spkLength,numres,p_val] = find_calciumdx_corr_pairs_2011(region,1,1000,0.05,'false',0.100);
if ~isfield(region,'userdata') || ~isfield(region.userdata,'corr')
    region.userdata.corr{1}.corr_pairs=corr_pairs;
    region.userdata.corr{1}.params.useGaussEvents=useGaussEvents;
    region.userdata.corr{1}.params.window_sec=win;
    region.userdata.corr{1}.params.spkLength_frames=spkLength;
    region.userdata.corr{1}.params.numres=numres;
    region.userdata.corr{1}.params.p_val=p_val;
else
    len=length(region.userdata.corr);
    region.userdata.corr{len+1}.corr_pairs=corr_pairs;
    region.userdata.corr{len+1}.params.useGaussEvents=useGaussEvents;
    region.userdata.corr{len+1}.params.window_msec=win;
    region.userdata.corr{len+1}.params.spkLength_frames=spkLength;
    region.userdata.corr{len+1}.params.numres=numres;
    region.userdata.corr{len+1}.params.p_val=p_val;
end
% useGaussEvents='false'; win=0.100; spkLength=5; numres=1000; p_val=0.05;
%data=region.userdata.corr{6}.corr_pairs{1};

data=corr_pairs{1};
fname2base = [fnm(1:end-4) 'dCorr_' num2str(win*1000) 'msWin-' num2str(p_val*100) 'perc' '-Gauss' useGaussEvents];
save([fname2base 'Region' '.mat' ],'region')

% dlmwrite([fname2base '.txt'],corr_pairs{1},'delimiter','\t','newline','unix');
%}

%{
d1 = data(1:10,1); d2 = data(1:10,2);
% adjMAT = zeros(length(region.contours));
% adjMAT(d1,d2) = 1;
% adjMAT = sparse(adjMAT);
% spy(adjMAT)

% DG = sparse([1 1 1 2 2 3 3 4 5 6 7 7 8 9 9  9 9],[2 6 8 3 1 4 2 5 4 7 6 4 9 8 10 5 3],true,10,10)

mx = max([d1; d2]);
DG = sparse(d1,d2,true,mx,mx)
h = view(biograph(DG));
[S,C] = graphconncomp(DG)
colors = jet(S);
for i = 1:numel(h.nodes)
  h.Nodes(i).Color = colors(C(i),:);
end
%}

%===Contour Graph plot of connectivity edges over image of cell contours as nodes===========
%myPlotImageCoords

%the following is important for getting the distances right if the data pixel dimensions are not equivalent
%And below the scripts will assume wherever 'rXY' is used, that it is szX (m dimension) which must be scaled up.
%the following assumes that the modulus of raster scanned data is 0 (equally divisible image size) and that for CCD images the ratio of image dimensions is either equivalent or not equally divisible
[szX,szY] = size(region.image);  %assuming szY is the largest dimension and szX may or may not need to be scaled.
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;
end

%the following sets which region.location you want to analyse
% numLoca = unique(region.location);
numLoca = 2;
figure;
set(gcf,'color',[1 1 1]);
% colormap(flipud(gray));
% colormap(ax(2),jet)
%plot gray dotted outline of each region----
hold on
h1=[];
for c = 1:size(region.contours,2)
    if region.location(c) == numLoca
    plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[0.5 0.5 0.5]);
%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0.7 0.7 0.7]);
%    set(h1(c),'FaceAlpha',0.1,'EdgeAlpha',0.1)
    else
    plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2)*rXY,'color',[1 1 1]);
%    h1(c)=patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
%    set(h1(c),'FaceAlpha',0,'EdgeAlpha',0)        
    end
end
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)]);
ylim([0 imagesize(1)*rXY]);
set(gca,'ydir','reverse');
set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
axis off

% i=4;
h2=[];
dist=[];
edgeList=[];
% i = 1;
distAll =[];
for i = 1:size(data,1)
   centr1 = centroid(region.contours{data(i,1)});
   centr2 = centroid(region.contours{data(i,2)});
   distAll = [distAll; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
end

for i = 1:size(data,1)
   centr1 = centroid(region.contours{data(i,1)});
   centr2 = centroid(region.contours{data(i,2)});
   if region.location(data(i,1)) == numLoca && region.location(data(i,2)) == numLoca
% %    h(i)=line([centr1(1) centr2(1)],[centr1(2) centr2(2)],'color','k');
%    h2(i)=patch([centr1(1) centr2(1)],[centr1(2) centr2(2)],[0 0 0]);
% %    h1(i) = patch(x,y,mycolors(numStim,:));
% %    set(h2(i),'EdgeColor',[1 0 0])
%    set(h2(i),'FaceAlpha',0.1,'EdgeAlpha',0.1)
   edgeList = [edgeList; data(i,:)];
   dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)*rXY-centr2(2)*rXY))^2)];
%    else
%        %    h(i)=line([centr1(1) centr2(1)],[centr1(2) centr2(2)],'color','k');
%    h2(i)=patch([centr1(1) centr2(1)],[centr1(2) centr2(2)],[1 1 1]);
% %    h1(i) = patch(x,y,mycolors(numStim,:));
%    set(h2(i),'FaceAlpha',0,'EdgeAlpha',0)
% %    dist = [dist; sqrt((abs(centr1(1)-centr2(1)))^2+(abs(centr1(2)-centr2(2)))^2)];
   end
end
% set(h2(i),'EdgeColor',[0 0 1])
% set(h1(i),'FaceAlpha',0.1,'EdgeAlpha',0.1)
% alpha(h,0.1);
mycolors = jet(numel(dist));
edgeList = [edgeList dist];
edgeList = sortrows(edgeList,3);
for i = 1:size(edgeList,1)
   centr1 = centroid(region.contours{edgeList(i,1)});
   centr2 = centroid(region.contours{edgeList(i,2)});
   h2(i)=patch([centr1(1) centr2(1)],[centr1(2)*rXY centr2(2)*rXY],[0 0 0]);
   set(h2(i),'EdgeColor',mycolors(i,:))
   set(h2(i),'FaceAlpha',0.5,'EdgeAlpha',0.5)
end
printfig('png',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.png'])
% plot2svg

% figure; hist(dist*region.spaceres,20,'color',[0.5 0.5 0.5])
% ylabel('No. of correlations'); xlabel('Distance (um)')
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',[0.5 0.5 0.5])






%======Normalized histogram of correlation distances=======================
%--get all pairwise distances between cell centroids in the image----
ct = zeros(length(region.contours),2);
for c = 1:size(ct,1)
    ct(c,:) = centroid(region.contours{c});
end

[a b] = meshgrid(ct(:,1)+sqrt(-1)*ct(:,2).*rXY);  %get all distances
dst = abs(a-b);
dst=dst.*region.spaceres;

figure; imagesc(dst); colorbar
printfig('epsc',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.eps'])

%--limit to unique pairwise distances ((N*N-1)/2)
dstLower = tril(dst,-1);
dstLower = dstLower(dstLower > 0);
% figure; hist(dstLower,20)
mean(dstLower)
mean(dist)
std(dstLower)/sqrt(numel(dstLower))
std(dist)/sqrt(numel(dist))
[h,p,ci]=ttest2(dist,dstLower)

%--set up the histogram
[y2,x2] = hist(dstLower,30);
[y1,x1] = hist(dist*region.spaceres,x2);
% bar(x1,y1,1)
% [y2, x2] = hist(dstLower,x1);

% bar(x2,y2,1)
y1Density = y1./sum(y1);  %relative frequency histogram (counts/(totalCounts)
y2Density = y2./sum(y2);  %relative frequency histogram (counts/(totalCounts)
binwidth = x1(2) - x1(1);
y1Density = y1Density./binwidth;  %relative density histogram (counts/(totalCounts*binwidth)
y2Density = y2Density./binwidth;  %relative density histogram (counts/(totalCounts*binwidth)
normDensity = y1Density./y2Density; %normalized frequency

figure; 
ax(1)=subplot(3,1,1);
% bar(x1,y1/trapz(x1,y1),1,'FaceColor',[0.5 0.5 0.5]);  %areal density
bar(x1,y1Density,1,'FaceColor',[0.5 0.5 0.5]); %relative frequency histogram
xlabel('Distance (um)'); ylabel('Density'); axis square
title('dist of pairwise cell correlations')
% xlim([0 300])
% ylim([0 0.1])
ax(2)=subplot(3,1,2);
% bar(x2,y2/trapz(x2,y2),1,'FaceColor',[0.5 0.5 0.5]); 
bar(x2,y2Density,1,'FaceColor',[0.5 0.5 0.5]); 
xlabel('Distance (um)'); ylabel('Density'); axis square
title('dst of all Pairwise cell distances')
% xlim([0 300])
% ylim([0 0.1])
ax(3)=subplot(3,1,3);
% bar(x2,y2/trapz(x2,y2),1,'FaceColor',[0.5 0.5 0.5]); 
bar(x2,normDensity,1,'FaceColor',[0.5 0.5 0.5]); 
xlabel('Distance (um)'); ylabel('Density'); axis square
title('dist/dst')
% xlim([0 300])
xlimits = xlim;
line(xlimits,[1 1],'Color','k','LineStyle','--')
linkaxes(ax(1:2),'y');
linkaxes(ax,'x');
zoom xon
printfig('epsc',[fname2base datestr(now,'yyyymmdd-HHMMSS') '.eps'])