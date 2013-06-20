st = get(dpdetectors,'value');
currdir = pwd;
cd(fullfile(calciumdxpath,'SignalDetectors'));
[ons, offs, param] = feval(mt(st).name(1:end-2),fnm,region);
cd(currdir)
region.onsets = ons;
region.offsets = offs;
region.detectorname = mt(st).name(1:end-2);
region.detectorparam = param;

set(bnext,'enable','on');

calciumdxPlotTrace;