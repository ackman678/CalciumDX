hold off
num = str2num(get(txcellnum,'string'));
if isempty(num)
    return
end
if num < 1
    num = size(nt,1);
    set(txcellnum,'string',num2str(num));
end
if num > size(nt,1)
    num = 1;
    set(txcellnum,'string',num2str(num));
end


set(cnt,'facecolor',[0 0 0]);
set(cnt(num),'facecolor',1-cl(region.location(num),:));

idx = region.transients(1,num);

set(radio(1),'value',0);
set(radio(2),'value',0);
set(radio(3),'value',0);
set(radio(4),'value',0);
set(radio(idx),'value',1);

trmenu = uicontextmenu;
uimenu(trmenu, 'Label', 'Add event', 'Callback', 'hevAddEvent');
plot(nt(num,:),'uicontextmenu',trmenu)

hold on
f = find(spk(num,:)==1);
g = find(dec(num,:)==1);
cmenu = zeros(1,length(f));
coffmenu = zeros(1,length(f));
for c = 1:length(f)
    cmenu(c) = uicontextmenu;
    uimenu(cmenu(c), 'Label', 'Delete event', 'Callback', ['selev = ' num2str(c) '; hevDeleteEvent;']);
    uimenu(cmenu(c), 'Label', 'Move onset', 'Callback', ['selev = ' num2str(c) '; hevMoveEvent;']);
    if c > 1
        uimenu(cmenu(c), 'Label', 'Combine with previous', 'Callback', ['selev = ' num2str(c) '; hevCombineEvent;']);
    end
    h = plot(f(c),nt(num,f(c)),'or','uicontextmenu',cmenu(c));
    
    coffmenu(c) = uicontextmenu;
    uimenu(coffmenu(c), 'Label', 'Move offset', 'Callback', ['selev = ' num2str(c) '; hevMoveEventOff;']);
    plot(g(c),nt(num,g(c)),'og','uicontextmenu',coffmenu(c));
end

xlim(xlimits);

set(gca,'buttondownfcn','hevZoom')
set(gcf,'KeyPressFcn','hevButtonDown')