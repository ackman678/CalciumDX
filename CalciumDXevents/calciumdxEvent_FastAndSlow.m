function [s d]=calciumdxEvent_FastAndSlow(tr,num,region)
% tr=tr;

% addedSign=size(region.traces,2); %addded by JA 13.06.07

% [sF dF] = calciumdxEvent_DetSingTrHannFast(tr,num,'no',region.onsets{num},region.offsets{num});
% [sS dS] = calciumdxEvent_DetSingTrHannSlow(tr,num,'no',region.onsets{num},region.offsets{num});
[sF dF] = calciumdxEvent_DetSingTrHannFast(tr,num,'no');
[sS dS] = calciumdxEvent_DetSingTrHannSlow(tr,num,'no');
[sU dU] = calciumdxdettrial(tr(num,:));
s=sF;
d=dF;
actFast=cell(0);
for i=1:length(sF)
    actFast{i}=(sF(i):dF(i));
end
actSlow=cell(0);
for i=1:length(sS)
    actSlow{i}=(sS(i):dS(i));
end
actUsua=cell(0);
for i=1:length(sU)
    actUsua{i}=(sU(i):dU(i));
end
% intersect fast and usual and use usual as onset
for i=1:length(sF)
    for j=1:length(sU)
        if length(intersect(actFast{i},actUsua{j}))>2
            s(i)=sU(j);
        end
    end
end
% % intersect fast and slow and use slow as offset
% for i=1:length(sF)
%     for j=1:length(sS)
%         if length(intersect(actFast{i},actSlow{j}))>2
%             d(i)=dS(j);
%         end
%     end
% end
% % intersect slow and usual and use usual as onset
% for i=1:length(sS)
%     for j=1:length(sU)
%         if length(intersect(actSlow{i},actUsua{j}))>2
%             sS(i)=sU(j);
%         end
%     end
% end
% add slow events if not detected
allFast=[];
for i=1:length(sF)
    allFast=[allFast actFast{i}];
end
for j=1:length(sS)
    if isempty(intersect(actSlow{j},allFast))
        s=[s sS(j)];
        d=[d dS(j)];
    end
end
% s

actAll=cell(0);
for i=1:length(s)
    actAll{i}=(s(i):d(i));
end
ss=[];
dd=[];
for i=1:length(s)-1
    for j=i+1:i+1
        if ~isempty(intersect(actAll{i},actAll{j}))
            uh=union(actAll{i},actAll{j});
            ss=[ss uh(1)];
            dd=[dd uh(end)];
        else
            ss=[ss s(i)];
            dd=[dd d(i)];
        end
    end
end
if length(s)>1
    if isempty(intersect(actAll{i},actAll{j}))
        ss=[ss s(j)];
        dd=[dd d(j)];
    end
end
ss=sort(ss);
dd=sort(dd);
s=[];
d=[];
i=1;
while i<length(ss)
    if ss(i)==ss(i+1)
        ss(i+1)=[];
        dd(i)=min([dd(i) dd(i+1)]);
        dd(i+1)=[];
    else
        i=i+1;
    end
end
i=1;
while i<length(ss)
    if dd(i)==dd(i+1)
        dd(i+1)=[];
        ss(i)=max([ss(i) ss(i+1)]);
        ss(i+1)=[];
    else
        i=i+1;
    end
end
s=ss;
d=dd;
i=1;
while i<=length(s)
    if d(i)-s(i)<2
        s(i)=[];
        d(i)=[];
    else
        i=i+1;
    end
end
i=1;
trac=(tr(num,:)-mean(tr(num,:)))/mean(tr(num,:))*100;
while i<=length(s)
    if abs(trac(d(i))-trac(s(i)))<3
        s(i)=[];
        d(i)=[];
    else
        i=i+1;
    end
end




%
% actFast=zeros(1,size(tr,2))
% for i=1:length(sF)
%     actFast(sF(i):dF(i))=1;
% end
% actSlow=zeros(1,size(tr,2))
% for i=1:length(sS)
%     actSlow(sS(i):dS(i))=1;
% end
% actUsua=zeros(1,size(tr,2))
% for i=1:length(sU)
%     actUsua(sU(i):dU(i))=1;
% end
%
% spk(num,s) = 1;
% dec(num,d) = 1;
% % figure(fig);
% % hold off;
% hevPlotTrace