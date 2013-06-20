%getSimpleWaveResults.m
filename=fnm;
% filename='/Users/ackman/Figures and images/2photon/100913/TSeries-09132010-1150-040_defaultwithFindPeaks.mat';

%number of ROIs
nrois=size(region.traces,1);

spk = zeros(size(region.traces));
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
end

%number of waves
nwaves=length(region.waveonsets);

%actv fraction
actv = sum(spk,2);
actvfraction = length(find(actv>0))/size(spk,1);

%wave actv fraction
waveactvfraction=mean(region.waveprops.waveactvfraction);

%wave frequency
wavefreq= length(region.waveonsets)/(size(region.traces,2).*region.timeres);

%wave isi
tmpres=region.timeres;
d=region.waveonsets;
e=region.waveoffsets;

if ~isempty(d)
d1 = [d size(region.traces,2)];
e1 = [0 e-d];
ints=diff([0 d1]) - e1;
if d(1) ~= 1
    ints=ints(2:end);
end
if d(end) ~= size(region.traces,2)
    ints=ints(1:end-1);
end
ints=ints*tmpres;
waveisi= mean(ints);
else
    waveisi=NaN;
end

%wave size
wavesizemicrons2=[];
for i=1:numel(region.waveprops.waveareapixels)
    wavesizemicrons2(i) = (sum(region.waveprops.waveareapixels{i}) * region.spaceres^2);
end
wavesize= mean(wavesizemicrons2(~isnan(wavesizemicrons2)));

%wave distance
wavedist=mean(region.waveprops.wavedistance_px(~isnan(region.waveprops.wavedistance_px)));

%wave speed
wavespeeds=[];
for wfr=1:length(region.waveonsets)
    consecutivespeeds=region.waveprops.wavespeeds_umpersec{wfr};
    wavespeeds(wfr)=mean(consecutivespeeds);
end
wavespeed=mean(wavespeeds(~isnan(wavespeeds)));

%wave direction
wavedir = mean(region.waveprops.wavedirection.theta_radians(~isnan(region.waveprops.wavedirection.theta_radians))) * (180/pi);


sprintf(['filename' '\t' 'nrois' '\t' 'nwaves' '\t' 'actvfraction' '\t' 'waveactvfraction' '\t' 'wavefreq.hz' '\t' 'm_waveisi.s' '\t' 'm_wavesize.um2' '\t' 'm_wavedist.um' '\t' 'm_wavespeed.umpersec' '\t' 'm_wavedir.degs'...
    '\r' filename '\t' num2str(nrois) '\t' num2str(nwaves) '\t' num2str(actvfraction) '\t' num2str(waveactvfraction) '\t' num2str(wavefreq) '\t' num2str(waveisi) '\t' num2str(wavesize) '\t' num2str(wavedist) '\t' num2str(wavespeed) '\t' num2str(wavedir)])

