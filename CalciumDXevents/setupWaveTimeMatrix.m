%=========================================================================================
function s = setupWaveTimeMatrix(region,param)
s = zeros(length(param.locationIndices),param.movieLength);
for i = 1:length(param.locationIndices)
	waveonsets = region.wavedata{param.locationIndices(i)}.waveonsets;
	s(i,waveonsets) = 1;
end
%s_thick = gauss_events(s,param.sig);
%figure; imagesc(s_thick); title('obs')
end
