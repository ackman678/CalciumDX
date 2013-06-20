function [tr, trhalo, param] = calciumdxTR_ReadTracesPrTIFFmedianfilt(fname,region,series1)
global t1 ff
% Program used by calciumdx
% Reads traces as average fluorescence values inside the contours
%This script is used for reading 
param = [];
% tfig = figure('Name','Reading traces...','NumberTitle','off','MenuBar','none','doublebuffer','on','units','normalized','position',[0.3    0.5    0.4    0.025]);

% inf = imfinfo(fname);
% 
% if strcmp(inf(1).ByteOrder,'little-endian');
%     readspec='ieee-le';
% elseif strcmp(inf(1).ByteOrder,'big-endian');
%     readspec='ieee-be';
% end
% fid = fopen(fname,'rb',readspec);

len=length(region.contours);
[szX,szY,szZ] = size(series1);

%implementation of median filter (default of 3x3) which is same as 1px radius median fitler in ImageJ.
d1=zeros([szX szY szZ],'uint16');
% d1=zeros([sz numPlanes]);
for i=1:szZ
   d1(:,:,i)=medfilt2(series1(:,:,i));  %use default medfilt2 from above (which has a 3x3 default). Median filter is the best one for getting rid of detector shot noise.
end
series1=d1;
clear d1



tr = zeros(len,szZ);
% prg = zeros(1,length(region.contours)+1);
% figure(tfig);
% subplot('position',[0 0 1 1]);
% set(gca,'xtick',[],'ytick',[]);

wtbar=waitbar(0,'Please wait...');
% figure;
% B=zeros(size(region.image));
for c = 1:len
        waitbar(c/len,wtbar);
%     disp(num2str(c))
    
    ps = round(region.contours{c});
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,region.contours{c}(:,1),region.contours{c}(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    f{c} = sub2ind([szX, szY],fy,fx);
    
 
%     B(f{c})=1;
%     imshow(B)
    
end
close(wtbar);

wtbar=waitbar(0,'Please wait...');
% prg = zeros(1,length(inf)+1);
for d = 1:szZ
    %     prg(d) = 1;
    %     figure(tfig);
    %     imagesc(prg);
    % d
    waitbar(d/szZ,wtbar);
    im=series1(:,:,d);
    %     drawnow
    for c = 1:length(region.contours)
        tr(c,d) = mean(im(f{c}));
    end
end
close(wtbar);


if region.halomode == 1
    trhalo = zeros(length(region.contours),size(series1, 1));
    halo = cell(1,length(region.contours));
    ff = [];
    for c = 1:length(region.contours)
        ff = [ff; f{c}];
    end
    ff = unique(ff);

    for c = 1:length(region.contours)
        ct = repmat(centroid(region.contours{c}),size(region.contours{c},1),1);
        halo{c} = (region.contours{c}-ct)*sqrt(1+region.haloarea)+ct;
        halo{c}(find(halo{c}(:,1)<1),1) = 1;
        halo{c}(find(halo{c}(:,2)<1),2) = 1;
        halo{c}(find(halo{c}(:,1)>size(region.image,2)),1) = size(region.image,2);
        halo{c}(find(halo{c}(:,2)>size(region.image,1)),2) = size(region.image,1);
    end

%     prg = zeros(1,length(region.contours)+1);
%     set(tfig,'Name','Reading halos...');
    for c = 1:length(region.contours)
%         prg(c) = 1;
%         figure(tfig);
%         imagesc(prg);
%         set(gca,'xtick',[],'ytick',[]);
%         drawnow

        ps = round(halo{c});
        [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
        inp = inpolygon(subx,suby,halo{c}(:,1),halo{c}(:,2));
        fx = subx(find(inp==1));
        fy = suby(find(inp==1));
        f{c} = setdiff(sub2ind(size(region.image'),fx,fy),ff);

        for d = 1:size(series1, 1)
            fseek(fid,inf(d).StripOffsets(1)+2*min(f{c})-2,'bof');
            [im count] = fread(fid,max(f{c})-min(f{c})+1,'*uint8');
            trhalo(c,d) = mean(im(f{c}-min(f{c})+1));
        end
    end
else
    trhalo = [];
end

% delete(tfig);
