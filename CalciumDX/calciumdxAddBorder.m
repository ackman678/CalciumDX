set(bord_add,'foreground',[1 0 0]);

stl = 1;
c = 1;
x = [];
y = [];
isb = 0;
[szM,szN] = size(region.image);

while stl
    [x(c) y(c) butt] = ginput(1);
    x(c) = round(x(c));
    y(c) = round(y(c));
    
	if x(c) < 1
		x(c) = 1;
	elseif x(c) > szN
		x(c) = szN;
	end
	if y(c) < 1
		y(c) = 1;
	elseif y(c) > szM
		y(c) = szM;
	end
    
    if butt > 1
        if c == 1
            isb = 1;
        else
            stl = 0;
        end
        if isb == 1
            [mn i] = min([x(c) maxx-x(c)+1 y(c) maxy-y(c)+1]);
            switch i
                case 1
                    x(c) = 1;
                case 2
                    x(c) = maxx;
                case 3
                    y(c) = 1;
                case 4
                    y(c) = maxy;
            end
        end
    end
    if c == 1
        h = plot(x,y,'yo');
    end
    if c == 2
        delete(h);
        h = plot(x,y,':+y');
    end
    if c > 2
        set(h,'xdata',x,'ydata',y);
    end
    c = c+1;
end
if isb == 0
    x = [x x(1)];
    y = [y y(1)];
    set(h,'xdata',x,'ydata',y);
end


bord{length(bord)+1} = [];
bord{end} = [get(h,'xdata')' get(h,'ydata')'];
        
bhand(end+1) = h;

set(bord_add,'foreground',[0 0 0]);

calciumdxDetermineRegions
