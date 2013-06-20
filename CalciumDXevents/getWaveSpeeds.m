function region = getWaveSpeeds(region)
%getWaveSpeeds
%James Ackman, 1/13/2011
%updated to work by location, 5/9/2011
locationMarkers = unique(region.location);
regionsAll = splitRegion(region);
for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    results = getWaveSpeedsByLocation(tmpregion,locationIndex);
    region.wavedata{locationIndex}.waveprops.wavespeeds_umpersec=results.wavespeeds_umpersec;
end

function results = getWaveSpeedsByLocation(region, locationIndex)
if isfield(region.wavedata{locationIndex},'wavecentroids')
    results.wavespeeds_umpersec={};
    
    for wfr=1:length(region.wavedata{locationIndex}.waveonsets)
        wavecentroid=region.wavedata{locationIndex}.wavecentroids{wfr};  %save data in region structure array
        disp(num2str(wfr)) %for testing
        %printout some results. If not enough detected waveframe centroids, no wavespeeds will be calculated.
        ind=1:size(wavecentroid,1);
        goodind=ind(find(~isnan(wavecentroid(:,1))));
        
        if length(goodind) > 1
%             diffind=diff(goodind)';  %iterative centroid frame times
            diffind=(goodind(2:end)-goodind(1))';  %speeds from 1st centroid frame
%             pxdist=sqrt((abs(diff(wavecentroid(goodind,1)))).^2 + (abs(diff(wavecentroid(goodind,2)))).^2);  %iterative moving centroid speeds
            pxdist=sqrt((wavecentroid(goodind(2:end),1)-wavecentroid(goodind(1),1)).^2 + (wavecentroid(goodind(2:end),2)-wavecentroid(goodind(1),2)).^2);  %speed from 1st centroid frame
            consecutivespeeds=(pxdist*region.spaceres)./(diffind*region.timeres);
            medianwavespeed=median(consecutivespeeds);  %median much more accurate than mean for representation of true wave speed here. Consecutive speeds usually not gaussian distributed, large outliers typically.
            meanwavespeed = mean(consecutivespeeds);
            wavespeed = meanwavespeed;  %this can be changed to median if for the data sets being analysed if necessary
            disp(['median wavespeed: ' num2str(medianwavespeed) ' um/sec'])   %for testing
            disp(['mean wavespeed: ' num2str(meanwavespeed) ' um/sec'])  %for testing
%             figure; hist(consecutivespeeds)  %for testing
%             disp(num2str(consecutivespeeds))  %for testing
            % pxdist=pxdist(~isnan(pxdist));
            wavepathlength = sum(pxdist) .* region.spaceres;
            disp(['total wavepathlength: ' num2str(wavepathlength) ' um'])  %for testing
            dxtime = (region.wavedata{locationIndex}.waveoffsets(wfr) - region.wavedata{locationIndex}.waveonsets(wfr)).*region.timeres;
            disp(['pathlength/waveduration: ' num2str(wavepathlength/dxtime) ' um/sec'])
            results.wavespeeds_umpersec{wfr}=consecutivespeeds;
        else
            disp('not enough detected center of masses')
            results.wavespeeds_umpersec{wfr}=NaN;
        end
        
    end
else
    disp('run getWaveCentroids.m first')
end