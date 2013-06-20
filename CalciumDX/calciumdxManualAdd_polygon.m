dummya = [];
butt = 1;
if get(shaperad(1),'value') == 1
    answer = inputdlg('ROI diameter (pixels):','Input for ROI',1,{'10'});
    r = str2double(answer{1});
    
    width=r; %width in pixels
    initialpos=[size(region.image,2)/2 size(region.image,1)/2];
    t = 0:pi/3:2*pi;
    xt=(sin(t) * width)+initialpos(1);
    yt=(cos(t) * width)+initialpos(2);
    h = impoly(imgax, [xt' yt']);
    fcn = makeConstrainToRectFcn('impoly',get(imgax,'XLim'),...
        get(imgax,'YLim'));
    setPositionConstraintFcn(h,fcn);
    % addNewPositionCallback(h,@(p) title(mat2str(p,3)));  %function callback example
    % h.addNewPositionCallback(@AppendPolyCoords);  %function callback example, but doesn't work
    % AppendPolyCoords(newcn,h.getPosition())
    % function newcn=AppendPolyCoords(newcn,cn)
    % newcn=[newcn cn];
    % end
    
    a='normal';
    while strcmp(a,'alt') == 0
        a=get(gcf,'SelectionType');
        if strcmp(a,'extend') == 1
            newcn = h.getPosition();
            
            delete(dummya)
            if size(newcn,1) < 3
                return
            end
            
            ct = calciumdxCentroid(newcn);
            reg = 0;
            for c = 1:length(region.coords)
                if inpolygon(ct(1),ct(2),region.coords{c}(:,1),region.coords{c}(:,2))
                    if reg == 0
                        reg = c;
                    elseif polyarea(region.coords{c}(:,1),region.coords{c}(:,2)) < polyarea(region.coords{reg}(:,1),region.coords{reg}(:,2))
                        reg = c;
                    end
                end
            end
            % lowar(reg) = tmp_str;
            
            if polyarea(newcn(:,1),newcn(:,2)) < lowar(reg) | polyarea(newcn(:,1),newcn(:,2)) > highar(reg)
                errordlg('Attempted contour area is not within limits!','Bad contour');
                return
            else
                cn{reg}{length(cn{reg})+1} = newcn;
            end
            
            tmp = num;
            num = reg;
            calciumdxDrawCells
            num = tmp;
        end
        waitforbuttonpress;
    end
    
    
    
else
    while butt <=1
        
        nx = [];
        ny = [];
        butt = 1;
        while butt <= 1
            [x y butt] = ginput(1);
            nx = [nx; x];
            ny = [ny; y];
            if size(nx,1) == 1
                dummya = plot(nx,ny,'+r');
            else
                set(dummya,'xdata',nx,'ydata',ny,'linestyle',':','linewidth',2,'marker','none');
            end
        end
        newcn = [nx ny];
        
        delete(dummya)
        if size(newcn,1) < 3
            return
        end
        
        ct = calciumdxCentroid(newcn);
        reg = 0;
        for c = 1:length(region.coords)
            if inpolygon(ct(1),ct(2),region.coords{c}(:,1),region.coords{c}(:,2))
                if reg == 0
                    reg = c;
                elseif polyarea(region.coords{c}(:,1),region.coords{c}(:,2)) < polyarea(region.coords{reg}(:,1),region.coords{reg}(:,2))
                    reg = c;
                end
            end
        end
        % lowar(reg) = tmp_str;
        
        if polyarea(newcn(:,1),newcn(:,2)) < lowar(reg) | polyarea(newcn(:,1),newcn(:,2)) > highar(reg)
            errordlg('Attempted contour area is not within limits!','Bad contour');
            return
        else
            cn{reg}{length(cn{reg})+1} = newcn;
        end
        
        tmp = num;
        num = reg;
        calciumdxDrawCells
        num = tmp;
    end
end
