dummya = [];
butt = 1;
% if get(shaperad(1),'value') == 1
%     answer = inputdlg('Cell radius (pixels):','Input for cell ROI',1,{'6'});
%     r = str2num(answer{1});
% end

% tmp_str = str2double(get(txarlow,'String'));

while butt <= 1
%     if get(shaperad(1),'value') == 1
%         x = [];
%         y = [];
%         [x(1) y(1) butt] = ginput(1);   %get first xy coord
%         if butt > 1
%             return
%         end
%         
%         %dummya = plot(x,y,'+r');
%         %  [x(2) y(2) butt] = ginput(1);  %get second xy coord
%         %    if butt > 1
%         %        delete(dummya)
%         %        return
%         %    end
%         
%         %r = sqrt((x(2)-x(1)).^2+(y(2)-y(1)).^2);  %pythagorean theorem to get radius
%         theta = 0:pi/50:2*pi-pi/50;
%         newcn = [r*cos(theta)'+x(1) r*sin(theta)'+y(1)];  %coords of the circle with the given radius r, at the point indicated by x1,y1
%     else
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
%     end
    
    delete(dummya)
    if size(newcn,1) < 3
        return
    end
    
    if allowRegionOverlap < 1
    
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
    
    else
        reg = num;
        
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