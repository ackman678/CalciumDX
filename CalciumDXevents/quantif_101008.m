tr = myfilter(region.traces(3,:),5);
figure; subplot(2,1,1)
plot(region.traces(3,:))
subplot(2,1,2)
plot(tr)
title('myfilter, hanning, order 5')
%----------------------------------------------------------------------------------------
wtbar=waitbar(0,'Please wait...');
for c = 1:len
    %     prg(c) = 1;
    %     figure(tfig);
    %     imagesc(prg);
    %     set(gca,'xtick',[],'ytick',[]);
    %     drawnow
    waitbar(c/len,wtbar);
%     im=series1{1,1};
    ps = round(region.contours{c});
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,region.contours{c}(:,1),region.contours{c}(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    f{c} = sub2ind([szX, szY],fy,fx);
end
close(wtbar);

tmp=zeros([256 512]);
figure;imagesc(tmp);
hold on;
plot(region.contours{43}(:,2),region.contours{43}(:,1),'.k')
plot(region.contours{43}(:,1),region.contours{43}(:,2),'.k')
%----------------------------------------------------------------------------------------

%get subindices row,col vectors for ROI no. 2
BWall=zeros(sz);
%         for i=1:num
i=2;
% [rows, cols] = find(L==i);
% %             region.contours{length(region.contours)+1} = [rows cols];
% %             region.contours{length(region.contours)+1} = [cols rows];
% region.contours{i}= [cols rows];

Atmp=zeros(sz);
%             Atmp(region.contours{i}(:,1),region.contours{i}(:,2))=1;
Atmp(region.contours{i}(:,2),region.contours{i}(:,1))=1;
figure; imshow(Atmp)

BWtmp=bwperim(Atmp);
BWall(BWtmp)=1;
figure; imshow(BWall)

figure;imshow(Atmp)
[Btmp] = bwboundaries(Atmp,'noholes');
hold on
plot(Btmp{1}(:,2), Btmp{1}(:,1), 'r', 'LineWidth', 2)

BWtmp=bwperim(Atmp);
[rows,cols]=ind2sub(sz,find(BWtmp))

BWall(BWtmp)=1;
%         end
figure; imshow(BWall)
        
%----------------------------------------------------------------------------------------
        
        figure;
%         subplot('position',[0.1 0.6 0.74 0.35]);
        imagesc(region.image)
        hold on;
        cl = hsv(length(region.name));
        cnt = zeros(1,length(region.contours));
        for c = 1:length(region.contours)
            cnt(c) = patch(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),[0 0 0]);
%             cnt(c) = patch(region.contours{c}(:,1),region.contours{c}(:,2),[1 1 1]);
            
%             cnt(c) = patch(region.contours{c}(:,2),region.contours{c}(:,1),[0 0 0]);
            set(cnt(c),'edgecolor',cl(region.location(c),:));
%             set(cnt(c),'ButtonDownFcn',['set(numslider,''value'',' num2str(c) '); calciumdxPlotTrace;']);
        end
%         axis equal
%         imagesize = size(region.image);
% %         xlim([0 imagesize(2)])
% %         ylim([0 imagesize(1)])
% %         set(gca,'ydir','reverse');
%         box on
%         set(gca,'color',[0 0 0]);
%         set(gca,'xtick',[],'ytick',[]);
%----------------------------------------------------------------------------------------

        figure;
        tmp=zeros(size(region.image));
        for i=1:length(region.contours) 
           imshow(tmp)
           plot(region.contours{i}(:,2),region.contours{i}(:,1))
           hold on
        end
        
        
        figure;
        BW = imread('blobs.png');
imshow(BW,[]);
s=size(BW);
for row = 2:55:s(1)
   for col=1:s(2)
      if BW(row,col),
         break;
      end
   end

   contour = bwtraceboundary(BW, [row, col], 'W', 8, 50,...
                                   'counterclockwise');
   if(~isempty(contour))
      hold on;
      plot(contour(:,2),contour(:,1),'g','LineWidth',2);
      hold on;
      plot(col, row,'gx','LineWidth',2);
   else
      hold on; plot(col, row,'rx','LineWidth',2);
   end
end

%----------------------------------------------------------------------------------------
myPlotRasterHist(fnm,region,[],'false','true')


%----------------------------------------------------------------------------------------
%hann filter
ntFilt=[];
for i=1:size(nt,1)
    y=myfilter(nt(i,:),2);
    ntFilt=[ntFilt; y];
end
figure; imagesc(ntFilt)

%moving average filter with window size 2 (movement artifacts result in zaxis changes involving 100% population dF change possible spanning 2 frames)
x=nt;
x=region.traces;
ntFilt2=[];
windowSize=2;
for i=1:size(x,1)
     y=filter(ones(1,windowSize)/windowSize,1,x(i,:)'); %same as rolling window mean filter of length 'windowSize'
    ntFilt2=[ntFilt2; y'];
end
figure; imagesc(ntFilt2)

%----------------------------------------------------------------------------------------

y=mean(ntFilt2,1); %mean of this moving average data
dy=diff(y);
dy(dy>0)=0;

figure; plot(dy);
dyAvg=mean(dy)
dySD=std(dy)
ySD=std(ntFilt2,1);

%can set offset threshold to 1 s.d of normal approx of median of data in calciumdxdettrial.m to get 
thr=median(abs(y))/0.6745; 
%or
thr=std(y);

%----------------------------------------------------------------------------------------
lambdahat = poissfit(data) 
figure();
plot(data(1,:));
hold;
plot(sig_pks1,data(1,sig_pks1),'or');
%----------------------------------------------------------------------------------------

%calciumdxevents help-- move artifact onset-------
fr=240; newfr=238;
[r,c]=find(region.artifactFrames == fr)
region.artifactFrames(r,c)=newfr;

%calciumdxevents help-- fuzzy delete all-----------------------------------------------------
region.artifactFrames=[76 90; region.artifactFrames]
region.artifactFrames=[790 810; region.artifactFrames]

badframes=[];
for i=1:size(region.artifactFrames,1)
    badframes{i}=region.artifactFrames(i,1):region.artifactFrames(i,2);
end

x=trSign*nt(num,:);
tr=x;
stn=s;
decpt=d;
    avgchange = mean(abs(diff(x)));
    stdchange = std(abs(diff(x)));

stnTemp=[1];
decptTemp=[1];
%         pknTemp=[1];
for c = 1:length(stn)
    stnTemp=[stnTemp stn(c)];
    decptTemp=[decptTemp decpt(c)];
    %             pknTemp=[pknTemp pkn(c)];
    for j=1:length(badframes)
        f=find(badframes{j}==stn(c), 1);
        if ~isempty(f)
            [dummy pk1]=min(tr(badframes{j}(end)+1:decpt(c)));
            pk1=pk1+badframes{j}(end);
            [dummy ind]= max(tr(badframes{j}(end)+1:pk1));
            ind=ind+badframes{j}(end);
            if ~isempty(pk1)
                if abs(tr(ind))+2*stdchange < abs(tr(pk1))
                    stnTemp(end)=ind;
                else
                    stnTemp=stnTemp(1:end-1);
                    decptTemp=decptTemp(1:end-1);
                    %                             pknTemp=pknTemp(1:end-1);
                end
            else
                stnTemp=stnTemp(1:end-1);
                decptTemp=decptTemp(1:end-1);
                %                         pknTemp=pknTemp(1:end-1);
            end
        end
    end
end

stn=stnTemp(2:end)
decpt=decptTemp(2:end)
%         pkn=pknTemp(2:end);
tr = region.traces;
%-----------------------------------------------------------------------------

%TSeries_038----------------------------------------------------
myPlotRasterHist(fnm,region,[],'false','true')
myPlotMultiTrace(region,[],spk,dec)
