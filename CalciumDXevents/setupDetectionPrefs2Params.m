function parameterArray = setupDetectionPrefs2Params(prefs)
parameterArray = {};
j=0;
for i = 1:length(prefs.params)
	j=j+1;
	parameterArray{j} = prefs.params(i).name;
	j=j+1;
	parameterArray{j} = prefs.params(i).value;
end
