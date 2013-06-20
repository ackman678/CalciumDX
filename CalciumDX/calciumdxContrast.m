brightness = get(bbright,'value');
contrast = get(bcontrast,'value');

contrast = contrast^0.25;
if contrast < 0.5
    m = 2*contrast/255;
else
    m = 1/(510-510*contrast+eps);
end
if brightness < 0.5
    b = -255*m*(1-(2*brightness)^0.25);
else
    b = 2*brightness-1;
end

cc = (0:255)*m+b;
cc(find(cc<0))=0;
cc(find(cc>1))=1;
climg = repmat(cc',1,3);

colormap(climg);