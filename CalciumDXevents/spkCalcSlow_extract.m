function [spkVet, ons, offs, minSpk]= spkCalcSlow_extract(data,time,deltaSp,mu,sigma,nStd);

Thresh=mu-nStd*sigma;
tempPrec=-inf;
ind= find(data< Thresh);
offs=[];
ons=[];
spkVet= [];
temp= [];
minSpk=[];
%          length(ind)
conk=0;
for k= ind
    %     conk=conk+1
    temp= time(k);
    if (temp - tempPrec > deltaSp)
        onsTemp=k;
        if onsTemp>1
            go='ok';
        else
            go='no';
        end
        if go=='ok'
        while data(onsTemp)<data(onsTemp-1) & go=='ok';
            onsTemp=onsTemp-1;
            if onsTemp<1
                go='no';
                onsTemp=1;
            end
        end
        end
        ons=[ons time(onsTemp)+1];
        offsTemp=k;
        if offsTemp<size(time,2)-2
            go='ok';
        else
            go='no';
        end
        if go=='ok'
        while data(offsTemp)<data(offsTemp+1) | data(offsTemp)<0 & go=='ok';
            offsTemp=offsTemp+1;
            if offsTemp>length(time)
                go='no';
                offsTemp=length(time);
            end
        end
        end
        offs=[offs time(offsTemp)];
        spkVet= [spkVet; temp];
        minSpk=[minSpk; data(k)];
    else
        if data(k)<minSpk(end)
            minSpk(end)=data(k);
            spkVet(end)=temp;
        end

    end
    tempPrec= temp;
end







