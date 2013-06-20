st = get(dpfilters,'value');
currdir = pwd;
cd(fullfile(calciumdxpath,'ImageFilters'));
[loca, param] = feval(mt(st).name(1:end-2),fnm,region);

region.filtername = mt(st).name(1:end-2);
region.filterparam = param;
cd(currdir)
num = 1;
calciumdxInputParams;