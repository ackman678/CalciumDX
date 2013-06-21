% function %hevArtifactToggle(hObject)
button_state = get(btoggle1,'Value');
if button_state == get(btoggle1,'Max')
 % Toggle button is pressed, take appropriate action
 showArtifacts=1;
 hevPlotTrace;
elseif button_state == get(btoggle1,'Min')
 % Toggle button is not pressed, take appropriate action
 showArtifacts=0;
 hevPlotTrace;
end