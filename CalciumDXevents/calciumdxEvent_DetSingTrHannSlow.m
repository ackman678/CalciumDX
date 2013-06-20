function [s, d] = calciumdxEvent_DetSingTrHannSlow(tr,nn,period,sss,ddd,region)



s=[];
d=[];
% [ss dd] = calciumdxdettrial(tr(nn,:)); %only used if lines 130-135 are used
% region.offsets{nn}=ddd; %only used for the plot lines below
% region.onsets{nn}=sss; %only used for the plot lines below

addedSign=1000;
lowFr=0.005;
freqSampl=1;
[b,a] = butter(2,2*lowFr/freqSampl,'low');
regSiz=size(tr,2);
regNum=[1:regSiz:size(tr,2)];
nstd=5;  % soglia x stimare rumore derivata
stdTh=3; % soglia su derivata
spkVet=[];

sign2an=[tr(nn,addedSign:-1:1) tr(nn,:) tr(nn,end:-1:end-addedSign+1)];
        sign2an = filter(b,a,sign2an);
% sign2an = myfilter(sign2an,5)-myfilter(sign2an,100);
%         trDif=tr(nn,:)-sign2an(addedSign+1:end-addedSign);
trDifAll=sign2an(addedSign+1:end-addedSign);
filSig=trDifAll;

% stdDevTm=[];
% for i=regNum
%     hhu=diff(trDifAll(i:i+regSiz-1),1);
%     avDer=mean(hhu);
%     stdDer=std(hhu);
%     ok=find(hhu>avDer-(1)*stdDer & hhu<avDer+(1)*stdDer);
%     stdDevTm=[stdDevTm std(hhu(ok))];
% %                     figure
% %                 [ui,uio]=hist(diff(trDifAll(i:i+regSiz-1),1),[-20:0.05:20]);
% %                 bar(uio,ui)
%
% end
%
% if max(stdDevTm)>min(stdDevTm)*1.5
%     regSiz=1000;
%     regNum=[1:regSiz:4000];
%     nstd=3;  % soglia x stimare rumore derivata
%     stdTh=5; % soglia su derivata
% else
%     regSiz=4000;
%     regNum=[1:regSiz:4000];
%     nstd=3;  % soglia x stimare rumore derivata
%     stdTh=1; % soglia su derivata
% end

for i=regNum
    addstd=-2;
    periodic='yes';
    while strcmp(periodic,'yes')
        periodic=period;
        addstd=addstd+2;%
        trDif=diff(trDifAll(i:i+regSiz-1));
        %         [trDif,r] = deconv(trDif,exp(-(1:30)/0.0));
        %         size(trDif)
        %         trDif=trDif+r(1:4000);

        %         trDif=smoothNeighbour(trDif,20);
        %         trDif=diff(trDif,2);
        % fit distribution
        numIntHist=200;
        %         if round((max(trDif)-min(trDif))/numIntHist)<20
        %             intHist=[min(trDif)-10:round((max(trDif)-min(trDif))/20):max(trDif)+10];
        %         else
        intHist=[min(trDif)-1:((max(trDif)-min(trDif))/numIntHist):max(trDif)+1];
        %         end
        [hcont,hx]=hist(trDif,intHist);
        intFit=hx(find(hcont>max(hcont)/4));
        avDer=mean(trDif(find(trDif>intFit(1) & trDif<intFit(end))));
        stdDer=std(trDif(find(trDif>intFit(1) & trDif<intFit(end))));
        sign2fit=trDif(  find(trDif>avDer-(nstd)*stdDer & trDif<avDer+(nstd)*stdDer) ); % sum(matAttOnOff)/numCell;
        [mu,sigma]=normfit(sign2fit);
        % end fit
% 
%                         figure
%                         [ui,uio]=hist(trDif,intHist);
%                         bar(uio,ui/(((max(trDif)-min(trDif))/numIntHist)*length(sign2fit)))
%                         hold on
%                         ny=normpdf(intHist,mu,sigma);
%                         plot(intHist,ny,'r-')
        
        
        [spkVetHelp, sh, dh, minSpk]= spkCalcSlow_extract(trDif,i:i+regSiz-2,3,mu,sigma,(stdTh+addstd));
        %       spkCalc_extract(trDif,i:i+regSiz-2,mu-(stdTh+addstd)*sigma, 3, -inf);
        %         [spkVetHelp, sh, dh, minSpk]= spkCalc_extract(trDif,i:i+regSiz-2,0, 3, -inf);
        s=[s sh];
        d=[d dh];
        %         [spkVetHelp2, tempPrec, minSpk]= spk_extract(-trDif,1:length(trDif),-mu-(stdTh+addstd)*sigma, 3, -inf);
        %         spkVetHelp=[spkVetHelp' spkVetHelp2'];
        %         spkVetHelp=sort(spkVetHelp);

        % periodic analysis
        [cont,x]=hist(diff(spkVetHelp),[0:5:30000]);

        [histSpk, tempPrec, minSpk]= spk_extract(-cont,1:length(cont),-(mean(cont(find(cont>0)))+1*std(cont(find(cont>0)))), 1, -inf);
        %         figure
        %         bar(x,cont)
        if isempty(histSpk)
            periodic='no';
        elseif histSpk(1)>15
            periodic='no';
        end
    end
    spkVet=[spkVet spkVetHelp' ];
end
spkVet=spkVet+1;
trDifAll=diff(trDifAll);
% s
% d
% s=[];
% d=[];
% for i=1:length(ss)
%     if isempty (find(spkVet>ss(i) & spkVet<dd(i)))==0
%         s=[s ss(i)];
%         d=[d dd(i)];
%     end
% end
% figure
% subplot(3,1,1)
% plot(tr(nn,:))
% set(gca,'xlim',[1 4000])
% subplot(3,1,2)
% % plot(filSig)
% % subplot(3,1,3)
% plot(trDifAll)
% set(gca,'xlim',[1 4000])
% hold on 
% plot([1 4000],mu-(stdTh+addstd)*sigma*ones(1,2))
% subplot(3,1,3)
% plot((0:size(tr,2)-1),100*(tr(nn,:)-mean(tr(nn,:)))/mean(tr(nn,:)));
% xlim([0 region.timeres*(size(tr,2)-1)])
% hold on
% plot((s-1),...
%     100*(tr(nn,s)-mean(tr(nn,:)))/mean(tr(nn,:)),'ro');
% plot((d-1),...
%     100*(tr(nn,d)-mean(tr(nn,:)))/mean(tr(nn,:)),'go');
% plot(spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-20,'k*')
% plot(spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-25,'r*')
% set(gca,'xlim',[1 4000])

%
% figure
% plot((trDifAll-mean(trDifAll))/max(trDifAll-mean(trDifAll))+3)
% hold on
% % plot(spkVet,ones(1,length(spkVet))*(mu-(stdTh+addstd)*sigma),'k*')
%
% plot((tr(nn,:)-mean(tr(nn,:)))/max(tr(nn,:)-mean(tr(nn,:))));
% plot((s-1),...region.offsets{nn}
%     (tr(nn,s)-mean(tr(nn,:)))/max(tr(nn,:)-mean(tr(nn,:)))-0.5,'ro');
% plot((d-1),...
%     (tr(nn,d)-mean(tr(nn,:)))/max(tr(nn,:)-mean(tr(nn,:)))-0.5,'go');
% plot(spkVet,(tr(nn,spkVet)-mean(tr(nn,:)))/max(tr(nn,:)-mean(tr(nn,:)))-0.5,'k*')
%
% plot((sign2an(addedSign+1:end-addedSign)-mean(sign2an(addedSign+1:end-addedSign)))/...
%     max(sign2an(addedSign+1:end-addedSign)-mean(sign2an(addedSign+1:end-addedSign)))-3)
%
% subplot(3,1,1)
% plot(trDifAll)
% hold on
% if isempty(spkVet)==0
%     plot(spkVet,mu-stdTh*sigma,'k*')
% end
% set(gca,'xlim',[5 4000])
% subplot(3,1,2)
% plot((0:size(tr,2)-1),100*(tr(nn,:)-mean(tr(nn,:)))/mean(tr(nn,:)));
% xlim([0 region.timeres*(size(tr,2)-1)])
% hold on
% plot((region.onsets{nn}-1),...region.offsets{nn}
%     100*(tr(nn,region.onsets{nn})-mean(tr(nn,:)))/mean(tr(nn,:)),'ro');
% plot((region.offsets{nn}-1),...
%     100*(tr(nn,region.offsets{nn})-mean(tr(nn,:)))/mean(tr(nn,:)),'go');
% plot(spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-20,'k*')
% set(gca,'xlim',[5 4000])
%
% subplot(3,1,3)
% plot(sign2an(addedSign+1:end-addedSign))
% figure
% plot((0:size(tr,2)-1),100*(tr(nn,:)-mean(tr(nn,:)))/mean(tr(nn,:)));
% xlim([0 region.timeres*(size(tr,2)-1)])
% hold on
% plot((s-1),...
%     100*(tr(nn,s)-mean(tr(nn,:)))/mean(tr(nn,:)),'ro');
% plot((d-1),...
%     100*(tr(nn,d)-mean(tr(nn,:)))/mean(tr(nn,:)),'go');
% plot(spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-20,'k*')
% plot(spkVet,100*(tr(nn,spkVet)-mean(tr(nn,:)))/mean(tr(nn,:))-25,'r*')
%
% set(gca,'xlim',[5 4000])
% % % figure
% %
% % %
