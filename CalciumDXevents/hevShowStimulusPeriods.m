% function %hevShowStimulusPeriods
button_state = get(btStimImport,'Value');
if button_state == get(btStimImport,'Max')
 % Toggle button is pressed, take appropriate action
 if isfield(region,'stimuli')
 showStimuli=1;
 hevPlotTrace;
 end
elseif button_state == get(btStimImport,'Min')
 % Toggle button is not pressed, take appropriate action
 showStimuli=0;
 hevPlotTrace;
end