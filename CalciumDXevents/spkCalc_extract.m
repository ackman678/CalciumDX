function [spkVet, ons, offs, minSpk]= spkCalc_extract(data,time,deltaSp,mu,sigma,nStd);

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
        go='ok';
        while data(onsTemp)<(mu-nStd*sigma) & go=='ok';
            onsTemp=onsTemp-1;
            if onsTemp<1
                go='no';
                onsTemp=1;
            end
        end
        ons=[ons time(onsTemp)+1];
        offsTemp=k;
        go='ok';
        while data(offsTemp)<(mu+1*sigma) & go=='ok';
            offsTemp=offsTemp+1;
            if offsTemp>length(time)
                go='no';
                offsTemp=length(time);
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








