delete(btProps1);
delete(btProps2);
delete(bord_add);
delete(bord_delete);
delete(bnext);

for c = 1:length(reg)
    region.coords{c} = reg{c};
    txlab(c) = uicontrol('Style','text','Units','normalized','String',['Region ' num2str(c)],'Position',[.87 .60-(c-1)*.07 .11 0.025],'FontSize',9,...
        'HorizontalAlignment','left','BackgroundColor',cl(c,:));
    inpt(c) = uicontrol('Style','edit','Units','normalized','String',['Name ' num2str(c)],'Position',[.87 .60-(c-1)*.07-0.035 .11 0.03],'FontSize',9,...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
end

bnext = uicontrol('Style','pushbutton','Units','normalized','String','Next >>','Position',[.87 .1 .05 .03],'FontSize',9, ...
    'Enable','on','Callback','calciumdxDetectCells');
