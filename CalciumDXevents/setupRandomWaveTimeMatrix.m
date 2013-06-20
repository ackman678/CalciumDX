%=========================================================================================
function s = setupRandomWaveTimeMatrix(param,s_obs)
%--set up random spike matrix dataset-------
s = zeros(size(s_obs));
for i = 1:param.nRegions
	if length(param.nRegions) < length(param.locationIndices) 
		s(i,:) = s_obs(i,:);  %set the first hemisphere to the obs hemsiphere wave times	
		i = i+1;   %if we are only reshuffling one hemisphere, then it will be the second ROI.
	end
	randWaveIntervals = getRandomWaveIntervals(param.nwaves(i), [], [], param.maxt,[],param.constrainNwaves);
	randWaveTimes=cumsum(randWaveIntervals);
	%disp(randWaveTimes)
	%disp(sum(randWaveIntervals))
	if ceil(randWaveTimes(end)/param.framePeriod) > param.movieLength && param.constrainNwaves > 0
		randWaveTimes = randWaveTimes-param.deltaspacing;
		s(i,ceil(randWaveTimes / param.framePeriod)) = 1;
	elseif param.constrainNwaves > 0
		s(i,ceil(randWaveTimes / param.framePeriod)) = 1;
	else
		wvIdx = find(randWaveTimes < param.maxt);
		s(i,ceil(randWaveTimes(wvIdx) / param.framePeriod)) = 1;
	end
	%error(num2str(randWaveTimes(end)))
end
%rand_sthick = gauss_events(s,param.sig);
%figure; imagesc(rand_sthick); title('res')
end