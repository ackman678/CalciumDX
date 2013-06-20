function [s, d] = calciumdxEvent_DetSingTrHPDer(tr,nn,period,region)



s=[];
d=[];
[ss dd] = calciumdxdettrial(tr(nn,:));

addedSign=size(tr,2);
lowFr=0.005;
freqSampl=1;
[b,a] = butter(2,2*lowFr/freqSampl,'high');
regSiz=10000;
regNum=[1:regSiz:4000];
nstd=3;  % soglia x stimare rumore derivata
stdTh=3; % soglia su derivata
spkVet=[];


for i=regNum
    addstd=-2;
    periodic='yes';
    while strcmp(periodic,'yes')
        periodic=period;
        addstd=addstd+2;%
        sign2an=[tr(nn,addedSign:-1:1) tr(nn,:) tr(nn,end:-1:end-addedSign+1)];
        filtSign = filter(b,a,sign2an);
        trDif=filtSign(addedSign+1:end-addedSign);
        trDif=diff(trDif);  %derivative
        % fit distribution
        numIntHist=200;
%         if round((max(trDif)-min(trDif))/numIntHist)<20
%             intHist=[min(trDif)-10:round((max(trDif)-min(trDif))/20):max(trDif)+10];
%         else
            intHist=[min(trDif)-10:round((max(trDif)-min(trDif))/numIntHist):max(trDif)+10];
%         end
        [hcont,hx]=hist(trDif,intHist);
        intFit=hx(find(hcont>max(hcont)/4));
        avDer=mean(trDif(find(trDif>intFit(1) & trDif<intFit(end))));
        stdDer=std(trDif(find(trDif>intFit(1) & trDif<intFit(end))));
        sign2fit=trDif(  find(trDif>avDer-(nstd)*stdDer & trDif<avDer+(nstd)*stdDer) ); % sum(matAttOnOff)/numCell;
        [mu,sigma]=normfit(sign2fit);
        % end fit

%         figure
%         [ui,uio]=hist(trDif,intHist);
%         bar(uio,ui/(round((max(trDif)-min(trDif))/numIntHist)*length(sign2fit)))
%         hold on
%         ny=normpdf(intHist,mu,sigma);
%         plot(intHist,ny,'r-')
% % 
        [spkVetHelp, tempPrec, minSpk]= spk_extract(trDif,1:length(trDif),mu-(stdTh+addstd)*sigma, 3, -inf);

        [cont,x]=hist(diff(spkVetHelp),[0:1:300]);

        [histSpk, tempPrec, minSpk]= spk_extract(-cont,1:length(cont),-(mean(cont(find(cont>0)))+1*std(cont(find(cont>0)))), 1, -inf);

        if isempty(histSpk)
            periodic='no';
        elseif histSpk(1)>15
            periodic='no';
        end
    end
    spkVet=[spkVet spkVetHelp'+i-1];
end
spkVet=spkVet+1;
trDifAll=trDif;


s=[];
d=[];
for i=1:length(ss)
    if isempty (find(spkVet>ss(i) & spkVet<dd(i)))==0
        s=[s ss(i)];
        d=[d dd(i)];
    end
end

% figure
% subplot(3,1,1)
% plot(trDifAll)
% hold on
% if isempty(spkVet)==0
%     plot(spkVet,mu-stdTh*sigma,'k*')
% end
% set(gca,'xlim',[5 4000])
% subplot(3,1,2)
% plot(region.timeres*(0:size(tr,2)-1),100*(tr(nn,:)-mean(tr(nn,:)))/mean(tr(nn,:)));
% xlim([0 region.timeres*(size(tr,2)-1)])
% hold on
% plot(region.timeres*(region.onsets{nn}-1),...
%     100*(tr(nn,region.onsets{nn})-mean(tr(nn,:)))/mean(tr(nn,:)),'ro');
% plot(region.timeres*(region.offsets{nn}-1),...
%     100*(tr(nn,region.offsets{nn})-mean(tr(nn,:)))/mean(tr(nn,:)),'go');
% plot(region.timeres*spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-20,'k*')
% set(gca,'xlim',[region.timeres*5 region.timeres*4000])
% 
% subplot(3,1,3)
% plot(region.timeres*(0:size(tr,2)-1),100*(tr(nn,:)-mean(tr(nn,:)))/mean(tr(nn,:)));
% xlim([0 region.timeres*(size(tr,2)-1)])
% hold on  
% plot(region.timeres*(s-1),...
%     100*(tr(nn,s)-mean(tr(nn,:)))/mean(tr(nn,:)),'ro');
% plot(region.timeres*(d-1),...
%     100*(tr(nn,d)-mean(tr(nn,:)))/mean(tr(nn,:)),'go');
% plot(region.timeres*spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-20,'k*')
% plot(region.timeres*spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-25,'r*')
% 
% set(gca,'xlim',[region.timeres*10 region.timeres*4000])
% figure

% 
