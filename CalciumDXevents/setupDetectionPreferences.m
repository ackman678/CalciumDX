function prefs = setupDetectionPreferences(region);

matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxprefs = fullfile(matlabUserPath,'calciumdxprefs.mat');
load(calciumdxprefs)

%--Make first preference set
prefs(1).name = 'normal';
prefs(1).params(1).name = 'hannfilterorder';
prefs(1).params(1).value = 2;
prefs(1).params(1).description = 'hannfilterorder';

prefs(1).params(2).name = 'sd';
prefs(1).params(2).value = 2;
prefs(1).params(2).description = 'Int: Number of std dev to make threshold for detection above baselineAverage';

prefs(1).params(3).name = 'sd2';
prefs(1).params(3).value = 1;
prefs(1).params(3).description = 'Int: no. of stddev above average local trace change';

prefs(1).params(4).name = 'sd3';
prefs(1).params(4).value = 1;
prefs(1).params(4).description = 'Int: No. of stddev above average whole trace change';

prefs(1).params(5).name = 'nonfilt';
prefs(1).params(5).value = 1;
prefs(1).params(5).description = 'logical: 0 | 1. Determines whether or not to use the orignal, non-low pass filtered data for signal onset refinement.';

prefs(1).params(6).name = 'hipass';
prefs(1).params(6).value = 'true';
prefs(1).params(6).description = 'str: "true" | "false". Use high pass filter for detrending baseline?';

prefs(1).params(7).name = 'slidingWinStartFrame';
prefs(1).params(7).value = 3;
prefs(1).params(7).description = 'Int: starting frame location for sliding window, def 3';

prefs(1).params(8).name = 'block_size';
prefs(1).params(8).value = 3;
prefs(1).params(8).description = 'Int: length of sliding window in frames in which to look for local peaks, def 3 frames';

prefs(1).params(9).name = 'start_baseline';
prefs(1).params(9).value = round(22/region.timeres);
prefs(1).params(9).description = 'Int: Baseline in frames for the sliding window from peak, def 22 sec';

prefs(1).params(10).name = 'end_baseline';
prefs(1).params(10).value = round(4.5/region.timeres);
prefs(1).params(10).description = 'Int: Baseline in frames for the sliding window from peak, def 4.5 sec';

prefs(1).params(11).name = 'blockSizeMultiplier';
prefs(1).params(11).value = 2;
prefs(1).params(11).description = 'Int: if next peak is inside N times the block size from last pk, then window will start from last peak';

prefs(1).params(12).name = 'windowAverage';
prefs(1).params(12).value = 'mean';
prefs(1).params(12).description = 'str: "mean" | "median". Use weighted mean or median for window baseline.';

prefs(1).params(13).name = 'baselineAverage';
prefs(1).params(13).value = 'all';
prefs(1).params(13).description = 'str: "window" | "all". F0 baseline, either the sliding window or whole trace';

prefs(1).params(14).name = 'maxOffsetTime';
prefs(1).params(14).value = round(5/region.timeres);
prefs(1).params(14).description = 'Int: Max no. frames for offset time, def 5 sec';


%--Make second preference set
prefs(2).name = 'waves';
prefs(2).params = prefs(1).params;
for i = 1:length(prefs(2).params)
	switch prefs(2).params(i).name
	 case 'hannfilterorder'
		 prefs(2).params(i).value = 10;
	 case 'nonfilt'
		 prefs(2).params(i).value = 0;
	 case 'hipass'
		 prefs(2).params(i).value = 'true';		
	 case 'maxOffsetTime'
		 prefs(2).params(i).value = 0;
	end
end


%--Add custom preference sets if they exist
nPrefs = length(prefs);
if exist('dxeventsPrefs','var')
	for i=1:length(dxeventsPrefs.detector)
		switch dxeventsPrefs.detector(i).name
		 case 'calciumdxdettrial'
			for j = 1:length(dxeventsPrefs.detector(i).prefs)
				nPrefs = nPrefs+1;
				prefs(nPrefs).name = dxeventsPrefs.detector(i).prefs(j).name;
				prefs(nPrefs).params = dxeventsPrefs.detector(i).prefs(j).params;
			end
		end
	end
end
