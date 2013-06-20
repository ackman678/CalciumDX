function [bigmesh] = makeMatchingROImeshgrid(region,locationIndex)
%2011-11 James B. Ackman
if nargin < 2 || isempty(locationIndex), locationIndex = [2 3]; end

if region.spaceres<3
    roiBin_width = 50;
elseif region.spaceres>3
    roiBin_width = 25;
end


bigmesh.contours = {};
bigmesh.location = [];
cn = cell(1,numel(region.coords));

minx = [];
for i = 1:length(locationIndex)
    dx(i)=abs(min(region.coords{locationIndex(i)}(:,1)) - max(region.coords{locationIndex(i)}(:,1)));
    minx = [minx min(region.coords{locationIndex(i)}(:,1))];
    dy(i)=abs(min(region.coords{locationIndex(i)}(:,2)) - max(region.coords{locationIndex(i)}(:,2)));
end

[mindx,mindxIdx] = min(dx);
[mindy,mindyIdx] = min(dy);
origin_x = min(minx);

[szX,szY] = size(region.image);  %assuming szY is the largest dimension
sz = [szX szY];
if mod(max([szY szX]),min([szY szX])) == 0
    rXY=szY/szX;
else
    rXY = 1;
end

strel_sz=[roiBin_width roiBin_width]; %new big ROI size in pixels. i.e. == 50pxx50px to give 113 µm x 113 µm for 2.27um/px

% locationIndex = [2 3];  %testing
Arr=zeros(sz);
for c = locationIndex
    if isempty(cn{c})
            cn{c}{1} = region.coords{c};
    end
    if ~isempty(cn{c})
        for d= 1:numel(cn{c})
            ps = round(cn{c}{d});
%             [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
%             inp = inpolygon(subx,suby,cn{c}{d}(:,1),cn{c}{d}(:,2));
%             fx = subx(inp==1);
%             fy = suby(inp==1);
%             f = sub2ind(size(region.image),fy,fx);
%             
%             B=zeros(size(region.image));
%             B(f)=1;
%                     figure; imshow(B)
%             
%             B2 = bwperim(B);
%             [I,J]=ind2sub(sz,find(B2==1));
%             %         figure; imshow(B2)
%             % hold on
            if region.orientation.value(1)< min(ps(:,2))
               shiftY = round(strel_sz(1)/2);
            else
                shiftY = 0;
            end
            
            I2 = (min(ps(:,2))+shiftY):strel_sz(1):(min(ps(:,2))+abs(diff(dy))+mindy);
%             J2 = origin_x:strel_sz(2):(minx+mindx);
            J2 = origin_x:strel_sz(2):(minx+abs(diff(dx))+mindx);


%             I2=min(I):strel_sz(1):max(I);
%             J2=min(J):strel_sz(2):max(J);
            
            I2arr=repmat(I2',1,numel(J2));
            J2arr=repmat(J2,numel(I2),1);
            
            Arr=zeros(sz); %<------------------------***
            Arr(I2arr,J2arr)=1;
%             figure; imshow(Arr)    %for TESTING   <---------------------
%             [i,j]=ind2sub(sz,find(Arr==1));
            
%             K=inpolygon(i,j,I,J);  %find those that are bounded by the original ROI image mask
            %         % plot(j(K),i(K),'.','color','white') %shows the pixels bounded by the mask.
            
%             A=zeros(sz);
%             idx=sub2ind(sz,i(K),j(K));
%             A(idx)=1;
            %         figure; imshow(A)
            %         % A(B2)=1;
            %         % figure; imshow(A)
            
            % SE = strel('octagon', 9);
            SE = strel('rectangle', [strel_sz(1)-1 strel_sz(2)-1]);

            Arr = imdilate(Arr,SE);
            figure; imshow(Arr)
            hold on;
            plot(cn{c}{d}(:,1),cn{c}{d}(:,2),'-b');
            hold off;
            
%             A2=imdilate(A,SE);
%             figure;imshow(A2)
            
            [L,num]=bwlabel(Arr);
            RGB = label2rgb(L);
%                     figure; imshow(RGB)  %for TESTING    <--------------
            
            
%                             figure; subplot(2,1,1)
%                             imshow(mat2gray(region.image))
%                             hold on
%                             plot(J,I,'.','color','red')
%                             subplot(2,1,2)
%                             imshow(RGB)
%                             hold on
%                             plot(J,I,'.','color','red')
            
            %get subindices row,col vectors for ROI no. 2
%             BWall=zeros(sz);
            for i=1:num
                %             region.location = [region.location d]; %numeric index for which region the contour is in.  Usually a vector of 1's the length of region contours
                bigmesh.location = [bigmesh.location c];
                [rows, cols] = find(L==i);
                %             region.contours{length(region.contours)+1} = [rows cols];
                
                Atmp=zeros(sz);
                Atmp(rows,cols)=1;
                %             Atmp(region.contours{i}(:,1),region.contours{i}(:,2))=1;
                %             Atmp(region.contours{i}(:,2),region.contours{i}(:,1))=1;
                %             BWtmp=bwperim(Atmp);
                %             [rows,cols]=ind2sub(sz,find(BWtmp));
                %             region.contours{length(region.contours)+1} = [cols rows];
                
                %             figure;imshow(Atmp)
                BWtmp = bwboundaries(Atmp,'noholes');
                bigmesh.contours{length(bigmesh.contours)+1} = [BWtmp{1}(:,2) BWtmp{1}(:,1)];
                %             hold on
                %             plot(BWtmp{1}(:,2), BWtmp{1}(:,1), 'r', 'LineWidth', 2)
                
                %             BWall(BWtmp)=1;
            end
            %         figure; imshow(BWall)
            %             end

        end
    else
        continue
    end
end
