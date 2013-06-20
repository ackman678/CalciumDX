% pause(.2)
str = sum(get(gcf,'currentcharacter'));

if str == 28 || str == 44
    num = num - 1;
    set(txcellnum,'string',num2str(num));
%     delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6); delete(h7); delete(h8); delete(h9); delete(h10);
%     clear trmenu cmenu coffmenu h1 h2 h3 h4 h5 h6 h7 h8 h9 h10
    hev2PlotTrace;
%     xlimits = [0 size(nt,2)+1];
end

if str == 29 || str == 46
    num = num + 1;
    set(txcellnum,'string',num2str(num));
%     delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6); delete(h7); delete(h8); delete(h9); delete(h10);
%     clear trmenu cmenu coffmenu h1 h2 h3 h4 h5 h6 h7 h8 h9 h10
    hev2PlotTrace;
%     xlimits = [0 size(nt,2)+1];
end
