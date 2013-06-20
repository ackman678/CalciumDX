function  num=dispCalcTracNetActComp(txcellnum, region, spk, dec, imgax, trax,nt)  % num
global bbright bcontrast
brightness = get(bbright,'value');
contrast = get(bcontrast,'value');

hold off
num = str2num(get(txcellnum,'string'));
if isempty(num)
    return
end
if num < 1
    num = size(nt,1);
    set(txcellnum,'string',num2str(num));
end
if num > size(nt,1)
    num = 1;
    set(txcellnum,'string',num2str(num));
end

f = find(spk(num,:)==1);
g = find(dec(num,:)==1);
region.onsets{num}=f;
region.offsets{num}=g;

% if numOld~=num
subplot(imgax)

hold off
imagesc(region.image)
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap gray
[maxy maxx] = size(region.image);
set(bbright,'value',brightness);
set(bcontrast,'value',contrast);
HippoContrast;
%     

cellSel=num;
for i=cellSel
    hold on;
    plot(region.contours{i}([1:end 1],1),region.contours{i}([1:end 1],2),'r-')
end
% end

subplot(trax)
hold off
plot(nt(num,:))
hold on

% hold on
[yLim]=get(gca,'ylim');
ymin=yLim(1);
ymax=yLim(2);
matAttOnOff=zeros(length(region.contours),length(region.traces));
matAttOn=zeros(length(region.contours),length(region.traces));
matAttOff=zeros(length(region.contours),length(region.traces));
numCell=length(region.contours);
conSpa=0;
for i=1:numCell
    ons=region.onsets{i};
    ofs=region.offsets{i};
    for j=1:length(ons)
        matAttOnOff(i,[ons(j):ofs(j)-1])=1;
        matAttOn(i,[ons(j)])=1;
        matAttOff(i,[ofs(j)])=1;
    end
end
unit=(ymax-ymin)/2;
offSetNet=ymin;
plot(sum(matAttOnOff)*unit/numCell+offSetNet,'k-')
plot(sum(matAttOn)*unit/numCell+offSetNet,'b-')
plot([1 size(region.traces,2)],ones(1,2)*(offSetNet+unit),'r:')
set(gca,'buttondownfcn','zoomDrugComp')

% xlim(xlimits);

% set(gcf,'KeyPressFcn','hevButtonDownNetAct')

% numOld=num;
