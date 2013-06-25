function [loca, param] = calciumdxIF_Localize(fname,region)

answer = inputdlg('Cell radius (pixels):','Input for the localization filter',1,{'5'});
if isempty(answer)
    loca = region.image;
    param = [];
    return
end
rad = str2num(answer{1});
param = rad;

tfig = figure('Name','Localizing...','NumberTitle','off','MenuBar','none','doublebuffer','on','units','normalized','position',[0.3    0.5    0.4    0.025]);

a = region.image;

%The following two commands pad the top and bottom of the image with extra rows %and columns from the sides and top and bottom of image equal to the number set %in rad (12 by default)
aat = [repmat(a(:,1),1,rad) a repmat(a(:,end),1,rad)];  
aat = [repmat(aat(1,:),rad,1); aat; repmat(aat(end,:),rad,1)];
%next we set up a blank matrix the size of region.image
loca = zeros(size(a,1),size(a,2));

pr = zeros(1,2*rad+2);
figure(tfig);
subplot('position',[0 0 1 1]);
set(gca,'xtick',[],'ytick',[]);
for c = -rad:rad
    pr(1,sum(pr)+1) = 1;
    figure(tfig);
    imagesc(pr);
    set(gca,'xtick',[],'ytick',[]);
    drawnow
    
    for d = -rad:rad
        loca = loca+aat(c+rad+1:end-rad+c,d+rad+1:end-rad+d);
        %figure(1);
        %imagesc(loca);
    end
end

loca = loca/(2*rad+1)^2;
loca = a./(loca+eps);

%see 'doc meshgrid' to understand the next two lines
[x y] = meshgrid(-rad:rad);
gs = exp(-(x.^2+y.^2)/rad);

%finally we perform cross-correlation of our filter over the image and remove %padding
loca = xcorr2(loca,gs);
loca = loca(rad+1:end-rad,rad+1:end-rad);

delete(tfig);