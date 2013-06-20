%Print drug comparison traces---------------------------------------------------------------------------
%to be used in conjunction with calciumdxDrugCompareTraces
%print trace for particular cell with subplot of high,low pass filter--------------------------------------
%Subplot of contour data-----------------------------------------------------
%JBA, Thursday, April 24, 2008 11:04 AM
figure();
%-----1
num=num1;
nt=nt1;
region=region1;
subplot(3,2,1);
hold on;
for c = 1:size(region.contours,2)
        plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
end
h = patch(region.contours{num}([1:end 1],1),region.contours{num}([1:end 1],2),'black');
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
%set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
title(fnm1(length(pathname)+1:end), 'Interpreter','none')

%------3
%figure();
xf=nt(num,:);
tres = region.timeres;
Nyq=0.5*(1/region.timeres);
%Wn = Df/(0.5 * Sf);  %Df = desired cutoff freq; Sf = sampling frequency (1/region.timeres)

%Adjust the Desired freq (Df) values in the following two lines for desired filter
%Filter order must be at least 3 times less than the signal length (e.g. for 1000 pt signal, filter order of 300 is ok)
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-2 may be needed for low pass filter. 10 is the default order of the low pass filter, but may need to be incr to 20 or higher for certain applications.
xf=filtfilt(fir1(300,0.01/Nyq,'high'),1,xf);  %need just for baseline correction typically
xf=filtfilt(fir1(10,1.5/Nyq,'low'),1,xf);  %1 or 1.5 works good

subplot(3,2,3);
plot((1:length(nt))*tres,xf,'color','k');
%xlabel('Time (s)'); 
ylabel('dF/F'); 
title(['cell ' num2str(num)]);

%------5
subplot(3,2,5);
plot((1:length(nt))*tres,nt(num,:),'color','k');
xlabel('Time (s)'); 
ylabel('dF/F'); 

%-----2
num=num1;
nt=nt2;
region=region2;
subplot(3,2,2);
hold on;
for c = 1:size(region.contours,2)
        plot(region.contours{c}([1:end 1],1),region.contours{c}([1:end 1],2),'color',[0.5 0.5 0.5]);
end
h = patch(region.contours{num}([1:end 1],1),region.contours{num}([1:end 1],2),'black');
axis equal
imagesize = size(region.image);
xlim([0 imagesize(2)])
ylim([0 imagesize(1)])
set(gca,'ydir','reverse');
%set(gca,'color',[0 0 0]);
set(gca,'xtick',[],'ytick',[]);
title(fnm2(length(pathname)+1:end), 'Interpreter','none')

%------4
%figure();
xf=nt(num,:);
tres = region.timeres;
Nyq=0.5*(1/region.timeres);
%Wn = Df/(0.5 * Sf);  %Df = desired cutoff freq; Sf = sampling frequency (1/region.timeres)

%Adjust the Desired freq (Df) values in the following two lines for desired filter
%Filter order must be at least 3 times less than the signal length (e.g. for 1000 pt signal, filter order of 300 is ok)
%usually 0.005 - 0.01 works for high pass filter. Anywhere from 0.1-2 may be needed for low pass filter. 10 is the default order of the low pass filter, but may need to be incr to 20 or higher for certain applications.
xf=filtfilt(fir1(300,0.01/Nyq,'high'),1,xf);  %need just for baseline correction typically
xf=filtfilt(fir1(10,1.5/Nyq,'low'),1,xf);  %1 or 1.5 works good

subplot(3,2,4);
plot((1:length(nt))*tres,xf,'color','k');
%xlabel('Time (s)'); 
ylabel('dF/F'); 
title(['cell ' num2str(num)]);

%------6
subplot(3,2,6);
plot((1:length(nt))*tres,nt(num,:),'color','k');
xlabel('Time (s)'); 
ylabel('dF/F'); 


figname = ['cell' num2str(num)]; fname2 = [fnm1(1:end-4) figname '.' 'eps']; saveas(gcf,fname2,'epsc');