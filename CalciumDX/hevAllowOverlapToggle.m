% function %hevAllowOverlapToggle(hObject)
button_state = get(btoggle1,'Value');
if button_state == get(btoggle1,'Max')
 % Toggle button is pressed, take appropriate action
 allowRegionOverlap=1;
elseif button_state == get(btoggle1,'Min')
 % Toggle button is not pressed, take appropriate action
 allowRegionOverlap=0;
end