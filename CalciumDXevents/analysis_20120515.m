%110208_04.mat
%First two frames deleted from tiff movie due to laser artifact-- must fix the stimulusParams saved in the file (shift to the left by two frame periods)
numStim = 1;
nstimuli = length(region.stimuli{numStim}.stimulusParams);
for i = 1:nstimuli
region.stimuli{numStim}.stimulusParams{i}.frame_indices = region.stimuli{numStim}.stimulusParams{i}.frame_indices - 2;
region.stimuli{numStim}.stimulusParams{i}.frame_times = region.stimuli{numStim}.stimulusParams{i}.frame_times - (2.*region.timeres.*1e06);
region.stimuli{numStim}.stimulusParams{i}.stimulus_times = region.stimuli{numStim}.stimulusParams{i}.stimulus_times - (2.*region.timeres.*1e06);
end
batchFetchStimResponseProps({filename},region,[],[-2000 5000]);


%TSeries-12202010-1424-002_motioncorrect2scanlinesFix2_medFilt-wavedet2.mat
%single cell analysis movie, that I had never finished putting through wavedetection
i = 1; locationIndex = 1; region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) = region.wavedata{locationIndex}.waveprops.wavedirection.theta_radians(i) - 3*(pi/4);