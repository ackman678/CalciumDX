set(bthide,'string','Hide');
ishid = 1;
numBKUP = num;

if isfield(region,'location')   %same info as isimported   %jba 2012-02-03
    locationMarkers = unique(region.location);
else
    locationMarkers = num;
end

for num = locationMarkers
    centr{num} = [];
    areas{num} = [];
    
    if isfield(region,'location')  %same t/f as isimported   %jba 2012-02-03
        fInd = find(region.location == num);
         for c = 1:length(fInd)
            centr{num}(c,:) = calciumdxCentroid(cn{num}{c});
            areas{num}(c) = polyarea(cn{num}{c}(:,1),cn{num}{c}(:,2));
        end
        
    else
        
        centr{num} = [];
        areas{num} = [];
        for c = 1:length(cn{num})
            centr{num}(c,:) = calciumdxCentroid(cn{num}{c});
            areas{num}(c) = polyarea(cn{num}{c}(:,1),cn{num}{c}(:,2));
        end
        
    end
    
    
    if ishandle(handl{num})
        delete(handl{num});
    end
    handl{num} = [];
    subplot('position',[0.02 0.02 0.82 0.96]);
    for c = 1:length(cn{num})
        h = plot(cn{num}{c}([1:end 1],1),cn{num}{c}([1:end 1],2),'color',cl(num,:),'LineWidth',1);
        handl{num} = [handl{num} h];
    end
    set(handl{num}(find(areas{num} < lowar(num) | areas{num} > highar(num))),'visible','off');
    
    
    
    %added coord outline 7/7/2011
    if ishandle(handlCoord{num})
        delete(handlCoord{num});
    end
    handlCoord{num} = [];
    hold on
    for numcoords = 1:length(region.coords)
        if prod(max(region.coords{numcoords})) ~= prod(size(region.image))
            hCoord = plot([region.coords{numcoords}(:,1); region.coords{numcoords}(1,1)], [region.coords{numcoords}(:,2); region.coords{numcoords}(1,2)],'--','color',[0.5 0.5 0.5]);
            handlCoord{num} = [handlCoord{num} hCoord];
        end
    end
end

f = findobj('Type','line','Visible','on');
if isempty(f)
    set(btdelete,'enable','off');
    set(btnextscr,'enable','off');
else
    set(btdelete,'enable','on');
    set(btnextscr,'enable','on');
end
zoom xon


contourplt = gca;
axes(contourplt)
num = numBKUP;