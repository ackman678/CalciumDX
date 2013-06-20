pilim(num) = str2num(get(txpilim,'string'));
for c = 1:length(handl{num})
    p = sum(sqrt(sum((cn{num}{c}([2:end 1],:)-cn{num}{c}(:,:)).^2,2)));
    ar = polyarea(cn{num}{c}(:,1),cn{num}{c}(:,2));
    api = p^2/(4*ar);
    if api > pilim(num)
        set(handl{num}(c),'linewidth',2);
    else
        set(handl{num}(c),'linewidth',1);
    end
end
refresh