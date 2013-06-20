st=get(dpreaders,'value');
if strcmp(mt(st).name(1:end-2),'calciumdxEvent_DetSingTrUsual')
    button='no';
else
    button = questdlg({'More strict detection?','More strict detection? '},'strict detection','yes','no','no');
end

prg = zeros(1,size(tr,1)+1);
tfigg = figure('Name','spike detection','NumberTitle','off','doublebuffer','on','units','normalized','position',[0.3    0.5    0.4    0.025]);
subplot('position',[0 0 1 1]);
set(gca,'xtick',[],'ytick',[]);
for c = 1:198%:size(tr,1);
    prg(c) = 1;
    figure(tfigg);
    imagesc(prg);
    set(gca,'xtick',[],'ytick',[]);
    drawnow
    [s d] = feval(mt(st).name(1:end-2),region,c,button);

%     [s d] = calciumdxdettrial(tr(c,:));
%     set(progtx,''String'',[''Detecting '' num2str(c) '' of '' num2str(size(nt,1))]);
    spk(c,:) = 0;
    dec(c,:) = 0;
    spk(c,s) = 1;
    dec(c,d) = 1;
end;
close (tfigg)
% set(progtx,''String'','''');
hevPlotTrace;
