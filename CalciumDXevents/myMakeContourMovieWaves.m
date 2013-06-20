%myMakeContourMovieWaves
function A = myMakeContourMovieWaves(fnm,region,includebaseline,includecolorbar)
%make new myMakeContourMovie script by outputting a multipage tiff that can be passed to ImageJ and mencoder for making a suitably sized avi file.
%James Ackman, 1/6/2011
%need your 'region' data structure loaded into workspace first--
% load region

%examples:
%A = myMakeContourMovieWaves([],region,'true','true')
%A = myMakeContourMovieWaves([],region)
%A = myMakeContourMovieWaves([pathname filename],region)

%START Make Movie----------------------------------------------------------
if nargin < 3, includebaseline='false'; end
if nargin < 4, includecolorbar='false'; end

% if isfield(region,'tracesFilt')
% region.traces=region.tracesFilt;
% end

[szX,szY] = size(region.image);  %assuming szY is the largest dimension
szZ = size(region.traces,2);
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
    szX=szY;  %to make the resulting images square, in case the data was raster scanned with less lines in one dimension--
else
    rXY = 1;
end

% figure; imagesc(A(:,:,1))
% hold on
% plot(region.contours{1}(:,1),region.contours{1}(:,2),'color','k')

%set up cell array the length of number of contours containing the pixels indices for region.image within each contour defined in region.contour using inpolygon
len=length(region.contours);
f={};
for c = 1:len
    %     prg(c) = 1;
    %     figure(tfig);
    %     imagesc(prg);
    %     set(gca,'xtick',[],'ytick',[]);
    %     drawnow
    %     waitbar(c/len,wtbar);
    %     im=series1{1,1};
    ps = round(region.contours{c});
    ps=[ps(:,1) rXY*ps(:,2)];
    if rXY > 1
        idx=find(ps(:,2) == min(ps(:,2)));
        if (min(ps(:,2))-rXY) > 0
            ps(idx,2)=min(ps(:,2))-rXY;
        else
            ps(idx,2) = 1;
        end
    end
    
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    %     inp = inpolygon(subx,suby,region.contours{c}(:,1),region.contours{c}(:,2));
    %     inp = inpolygon(subx,suby,region.contours{c}(:,1),rXY*region.contours{c}(:,2));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    %     disp(num2str(c))
    f{c} = sub2ind([szX, szY],fy,fx);
end

%---------------------------------------------
%--Build script to intersect nt trace image values with events-------------
%----must first do hipass filtering for baseline correction of region.traces to improve the intersected nt trace image of events (otherwise the colormap is scaled to the first waves, if there are strong baseline changes in mean fluorescence)---
Nyq=0.5*(1/region.timeres);
hipasscutoff = 0.005;
ntFilt=zeros(size(region.traces));

if 900 > size(region.traces,2)
    highfilterorder=round((size(region.traces,2)/3)) - 1;  %data must of length of greater than 3x filter order.
elseif size(region.traces,2) == 900
    highfilterorder=round((size(region.traces,2)/3)) - 2;
else
    highfilterorder = 300;
end

for i=1:size(region.traces,1)
    %     xf=filtfilt(fir1(300,hipasscutoff,'high'),1,region.traces(i,:));
    xf=filtfilt(fir1(highfilterorder,hipasscutoff,'high'),1,region.traces(i,:));
    ntFilt(i,:)=xf;
end
ntFilt=mat2gray(ntFilt);

%------Now intersect the new ntFilt values with event times for a thresholded event matrix coded by event amplitudes
% figure; imshow(mat2gray(ntFilt)); colormap(jet)
ntSpk=zeros(size(region.traces));
for c=1:size(region.traces,1)
    if ~isempty(region.onsets{c})
        for d=1:length(region.onsets{c})
            %            ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = region.traces(c,region.onsets{c}(d):region.offsets{c}(d));
                        ntSpk(c,region.onsets{c}(d):region.offsets{c}(d)) = ntFilt(c,region.onsets{c}(d):region.offsets{c}(d));
        end
    end
end

% figure; imagesc(ntSpk)
% figure; imshow(mat2gray(ntSpk)); colormap(flipud(gray))   %this flips the colormap so we have a white background
% figure; imshow(mat2gray(ntSpk)); colormap(jet)
%---------------------------------------------
%{
s = zeros(len,szZ);
for c = 1:len
   for d = 1:length(region.onsets{c})
      s(c,region.onsets{c}(d):region.offsets{c}(d)) = 1;
   end
end

tr = region.traces;
nt = [];
for c = 1:size(tr,1)
    nt(c,:) = dfoverf(tr(c,:))*100;
end
% tr=nt;
%}

%{
A=zeros([szX szY szZ]);
for i=1:szZ
    Adummy=zeros(szX,szY);
    for c=1:len
        if s(c,i) == 1
        Adummy(f{c})=tr(c,i);
        end
    end
    A(:,:,i)=Adummy;
end
%}
m1 = mean(ntSpk,1);

if strcmp(includebaseline,'true')
    %     A=zeros([szX szY szZ]);
    %     for i=1:szZ
    %         Adummy=zeros(szX,szY);
    %         for c=1:len
    %             Adummy(f{c})=ntFilt(c,i);  %set ntSpk to ntFilt if you want a contour movie that includes baseline fluctuations in between detected events.
    %             if m1(i) > 0 && strcmp(includecolorbar,'true')
    %                 idx=find(ntSpk(:,i) > 0);
    %                 colorbarvalue=mean(ntSpk(idx,i));
    %                 Adummy(2:40,2:10)=colorbarvalue;
    %             end
    %         end
    %         A(:,:,i)=Adummy;
    %     end
    A=zeros([szX szY szZ],'uint8');
    ntSpk8bt = im2uint8(ntSpk);
    ntFilt8bt = im2uint8(ntFilt);
    for i=1:szZ
        Adummy=zeros(szX,szY,'uint8');
        for c=1:len
            Adummy(f{c})=ntFilt8bt(c,i);  %set ntSpk to ntFilt if you want a contour movie that includes baseline fluctuations in between detected events.
            if m1(i) > 0 && strcmp(includecolorbar,'true')
                idx=find(ntSpk8bt(:,i) > 0);
                colorbarvalue=mean(ntSpk8bt(idx,i));
                Adummy(2:40,2:10)=colorbarvalue;
            end
        end
        A(:,:,i)=Adummy;
    end
else
    %     A=zeros([szX szY szZ]);
    %     for i=1:szZ
    %         Adummy=zeros(szX,szY);
    %         for c=1:len
    %             Adummy(f{c})=ntSpk(c,i);  %set ntSpk to ntFilt if you want a contour movie that includes baseline fluctuations in between detected events.
    %         end
    %         A(:,:,i)=Adummy;
    %     end
    
    A=zeros([szX szY szZ],'uint8');
    ntSpk8bt = im2uint8(ntSpk);
    for i=1:szZ
        Adummy=zeros(szX,szY,'uint8');
        for c=1:len
            Adummy(f{c})=ntSpk8bt(c,i);  %set ntSpk to ntFilt if you want a contour movie that includes baseline fluctuations in between detected events.
        end
        A(:,:,i)=Adummy;
    end
    
end

% A=mat2gray(A);
% implay(A)

%{
%if you want to add contours for the cells, use the following code-- pretty slow because of bwperim inside the second for loop (much faster if only inside first loop, but then it can't distinguish between neighboring cells which have some overlapping pixels).
for i=1:szZ
    Adummy=zeros(szX,szY);
    B2=zeros(szX,szY);
    for c=1:len
        Adummy(f{c})=1;
        B = bwperim(Adummy);
        B2(B)=1;
    end
    %     B = bwperim(Adummy);
    %     B=mat2gray(B);
    Atemp=A(:,:,i);
    Atemp=Atemp+B2;
    A(:,:,i)=Atemp;
end
implay(A)
% implay(mat2gray(A))

Adummy=zeros(szX,szY);
Adummy(f{1})=1;
outline=sub2ind([szX szY],round(region.contours{1}(:,2)),round(region.contours{1}(:,1)));
Adummy(outline)=1;
figure; imshow(mat2gray(Adummy))

%use bwperim to find perim of ROI
B = bwperim(Adummy);
figure; imshow(mat2gray(B))

%use convhullto find perim of ROI instead?  No... get pixel indices for convex hull, but doesn't necessarily include all perimeter pixels
Adummy=zeros(szX,szY);
[I,J]=ind2sub([szX szY],f{1});
k=convhull(I,J);
outline=sub2ind([szX szY],I(k),J(k));
Adummy(outline)=1;
figure; imshow(mat2gray(Adummy))

%}

% movie2avi(mov,'quantif_100925_motioncorrectAvgtemplate8scanlines.avi')
%END Make Movie

if ~isempty(fnm)
    %Write to multipage TIFF
    fnm2=fnm(1:end-4);
    fnm2=[pathname fnm2];
    if strcmp(includebaseline,'true')
        fnm2 = [fnm2 '_contourmovie2withbackground.tif'];
    else
        fnm2 = [fnm2 '_contourmovie2.tif'];
    end
    % fnm2=[pathname 'TSeries_038_multipage.tif'];
    % fnm2=['TSeries_038_multipage.tif'];
    % fnm2=['TSeries_040_multipage.tif'];
    for i = 1:size(A,3)
        imwrite(A(:,:,i), fnm2,'tif', 'Compression', 'lzw', 'WriteMode', 'append');
        %         imwrite(A(:,:,i), fnm2,'tif', 'Compression', 'none', 'WriteMode', 'append');
    end
end