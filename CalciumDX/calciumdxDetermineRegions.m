if length(bord) > 0
    set(bord_delete,'enable','on');
    reg = cell(1,0);
	reg{end+1} = [1 1; 1 maxy; maxx maxy; maxx 1];
	for c = 1:length(bord)
        reg{end+1} = bord{c}(1:end-1,:);
    end
else
    reg = cell(1,1);
    reg{1} = [1 1; 1 maxy; maxx maxy; maxx 1];
    set(bord_delete,'enable','off');
end

%{
if length(bord) > 0
    set(bord_delete,'enable','on');
    pt = [];
    for c = 1:length(bord)
        pt = [pt; [bord{c}(1,:) c]];
        pt = [pt; [bord{c}(end,:) c]];
    end
    lst = [];
    
    lst = [1 1 0];
    f = find(pt(:,1)==1);
    lst = [lst; sortrows(pt(f,:))];
    lst = [lst; [1 maxy 0]];
    f = find(pt(:,2)==maxy);
    lst = [lst; sortrows(pt(f,:))];
    lst = [lst; [maxx maxy 0]];
    f = find(pt(:,1)==maxx);
    lst = [lst; flipud(sortrows(pt(f,:)))];
    lst = [lst; [maxx 1 0]];
    f = find(pt(:,2)==1);
    lst = [lst; flipud(sortrows(pt(f,:)))];
    
    usd = zeros(1,size(lst,1));
    reg = cell(1,0);
    
    if size(lst,1) == 4
        reg{end+1} = [1 1; 1 maxy; maxx maxy; maxx 1];
    end
    
    for c = 1:size(lst,1)
        if usd(c) == 0 & lst(c,3) > 0
            reg{end+1} = zeros(0,3);
            j = c;
            while size(reg{end},1) < 2 | reg{end}(1,3) ~= reg{end}(end,3)
                usd(j) = 1;
                reg{end} = [reg{end}; lst(j,:)];
                j = mod(j,size(lst,1))+1;
                if lst(j,3) > 0
                    reg{end} = [reg{end}; lst(j,:)];
                    if bord{lst(j,3)}(1,1) == reg{end}(end,1) & bord{lst(j,3)}(1,2) == reg{end}(end,2)
                        reg{end} = [reg{end}; [bord{lst(j,3)}(2:end-1,:) lst(j,3)*ones(size(bord{lst(j,3)},1)-2,1)]];
                    else
                        reg{end} = [reg{end}; flipud([bord{lst(j,3)}(2:end-1,:) lst(j,3)*ones(size(bord{lst(j,3)},1)-2,1)])];
                    end
                    j = setdiff(find(lst(:,3)==lst(j,3)),j);
                end
            end
            reg{end} = reg{end}(:,1:2);
        end
    end
    
    for c = setdiff(1:length(bord),lst(:,3))
        reg{end+1} = bord{c}(1:end-1,:);
    end
    
else
    reg = cell(1,1);
    reg{1} = [1 1; 1 maxy; maxx maxy; maxx 1];
    set(bord_delete,'enable','off');
end
%}

regax = subplot('position',[0.87 0.64 0.11 0.10]);
delete(get(regax,'children'));
hold on
cl = hsv(length(reg));
for c = 1:length(reg)
    patch(reg{c}(:,1),reg{c}(:,2),cl(c,:));
end
set(gca,'xtick',[],'ytick',[],'ydir','reverse');
axis equal
axis tight
box on
