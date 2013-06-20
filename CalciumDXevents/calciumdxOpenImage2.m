[filename, pathname] = uigetfile({'*.tif'}, 'Choose image to open');
fnm = [pathname filename];

delete(logoim)
delete(logoname)
delete(myname)

imgax = subplot('position',[0.02 0.02 0.82 0.96]);
imagesc(a);
hold on

set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on

colormap gray

set(bzoom,'enable','on');
set(bbright,'enable','on');
set(bcontrast,'enable','on');
set(bnext,'enable','on');
set(inptsr,'enable','on');
set(inpttr,'enable','on');


[maxy maxx] = size(a);

zoom on;
calciumdxContrast;