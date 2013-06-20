wd = zeros(1,length(handl{num}));
for c = 1:length(handl{num})
    wd(c)= get(handl{num}(c),'linewidth');
end
f = find(wd==2);

tempcn = [];
spl = [];
for c = 1:length(cn{num})
    if isempty(find(f==c))
        tempcn{length(tempcn)+1} = cn{num}{c};
    else
        newcn = {};    
        crd = cn{num}{c};
        x = max([fix(min(crd(:,1)))-2 1]):min([fix(max(crd(:,1)))+2 size(a,2)]);
        y = max([fix(min(crd(:,2)))-2 1]):min([fix(max(crd(:,2)))+2 size(a,1)]);
        vls = loca(y,x);
        [xs ys] = meshgrid(x,y);
        in = inpolygon(xs,ys,crd(:,1),crd(:,2));
        vls = vls.*in;
        
        mx = zeros(size(vls));
        for xf = -1:1
            for yf = -1:1
                mx(2:end-1,2:end-1) = max(cat(3,mx(2:end-1,2:end-1),vls((2:end-1)+yf,(2:end-1)+xf)),[],3);
            end
        end
        
        [j i] = find(vls>=mx & vls~=0);
        i = x(1)+i-1;
        j = y(1)+j-1;
        
        dst = [];
        for d = 1:length(i)
            dst(:,d) = sum((crd-repmat([i(d) j(d)],size(crd,1),1)).^2,2);        
        end
        [mn bestcell] = min(dst,[],2);
        set(handl{num}(c),'visible','off');
        for d = 1:length(i)
            newcn{d} = crd(find(bestcell==d),:);
            v1 = newcn{d}([2:end 1],:)-newcn{d};
            v2 = newcn{d}([end 1:end-1],:)-newcn{d};
            angl = sum(v1.*v2,2)./(sum(v1.^2,2).*sum(v2.^2,2)+eps);
            newcn{d} = newcn{d}(find(angl<0),:);
            if ~isempty(newcn{d})
                spl = [spl plot(newcn{d}([1:end 1],1),newcn{d}([1:end 1],2),'linewidth',2,'Color',1-cl(num,:))];
            end
            tempcn{length(tempcn)+1} = newcn{d};
        end
        drawnow
        refresh
    end
end

cn{num} = [];
for c = 1:length(tempcn)
    if polyarea(tempcn{c}(:,1),tempcn{c}(:,2)) > 0
        cn{num}{size(cn{num},2)+1} = tempcn{c};
    end
end

isadjust(num) = 1;

delete(spl)
calciumdxDrawCells