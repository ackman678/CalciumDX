function strel_sz = getStrelsize(region,c)
%get structured element size dimensions from defined ROIs
%--the following gives us the square dimensions in pixels of the ROIs (strel_sz)-----
    if nargin < 2; c=1; end
    f={};
    %the following assumes that the modulus of raster scanned data is 0 (equally divisible image size) and that for CCD images the ratio of image dimensions is either equivalent or not equally divisible
    [szX,szY] = size(region.image);  %assuming szY is the largest dimension and szX may or may not need to be scaled.
%     szZ = size(region.traces,2);
    if mod(max([szY szX]),min([szY szX])) == 0
        rXY=szY/szX;
        szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
    else
        rXY = 1;
    end
    
    ps = round(region.contours{c});
    ps=[ps(:,1) rXY*ps(:,2)];
    if rXY > 1
        idx=find(ps(:,2) == min(ps(:,2)));
        ps(idx,2)=min(ps(:,2))-rXY;
    end
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    %f{c} = sub2ind([szX, szY],fy,fx);
    strel_sz=[numel(unique(fx)) numel(unique(fx))];