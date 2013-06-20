figure();
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP_HF(trSign*region.traces,num,'no',trNew); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;


figure();
ax(1)=subplot(2,1,1)
% tr = myfilter(region.traces(num,:),2);
dx=diff(region.traces(num,:));
% plot(tr);
plot([dx mean(dx)]);
ax(2)=subplot(2,1,2)
% spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxEvent_DetSingTrHP_HF(trSign*region.traces,num,'no',trNew); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;
spk(num,:) = 0; dec(num,:) = 0; [s d] = calciumdxdettrial(trSign*region.traces(num,:)); spk(num,s) = 1; dec(num,d) = 1; if region.transients(1,num) == 1 && sum(spk(num,:)) > 0; region.transients(1,num) = 4; end; hev2PlotTrace;
linkaxes(ax,'x');


data = clipboard('pastespecial');
roi1=(data.A_pastespecial(:,2))';

bkgrd=(data.A_pastespecial(:,2))';
roi1=region.traces(num,:);


tmp1=roi1;
figure(); 
ax(1)=subplot(3,1,1);
plot(tmp1)
ax(2)=subplot(3,1,2);
plot(tmp1-bkgrd)
ax(3)=subplot(3,1,3);
plot(bkgrd)
linkaxes(ax,'x');
set(gca,'UserData',[0 size(tmp1,2)])
set(gca,'buttondownfcn','hevZoom2')


profile on
calciumdxevents3
profile viewer
p = profile('info');
profsave(p,'profile_results')
