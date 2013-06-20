function [spkVet, tempPrec, minSpk]= spk_extract(data,time,Thresh, deltaSp, tempPrec);

ind= find(data< Thresh);
	    spkVet= [];
         temp= [];
         minSpk=[];
%          length(ind)
         conk=0;
for k= ind    
%     conk=conk+1
   temp= time(k);
       if (temp - tempPrec > deltaSp)
            spkVet= [spkVet; temp];
            minSpk=[minSpk; data(k)]; 
       elseif data(k)<minSpk(end)
            minSpk(end)=data(k);
            spkVet(end)=temp;
       end
   tempPrec= temp;
end
         
        
         
         

         

