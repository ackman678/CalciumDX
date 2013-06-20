function allH=showAllCont(region1,imgax1)
allH=[];
subplot(imgax1);
for i=1:length(region1.contours);
    hold on;
    h=plot(region1.contours{i}([1:end 1],1),region1.contours{i}([1:end 1],2),'r-');
    allH=[allH h];
end