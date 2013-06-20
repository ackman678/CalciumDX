function [A2, A3] = dilateContours(A,strel_sz)
%function to dilate structured elements (ROIs) in an array
%2011-11 James B. Ackman
A2 = zeros(size(A),'uint8');
A3 = zeros(size(A),'uint8');
hbar = waitbar(0,'Please wait...');
szZ = size(A,3);
for fr = 1:szZ
    
    Afr = A(:,:,fr);
    
    if ~isempty(find(Afr(:), 1))
        
        SE = strel('rectangle', [2 2]);
        Afr=imdilate(Afr,SE);
        %use rectangle strel object to dilate image to fill in gaps of 1 ROI in size.  e.g. 22x22 for 19x19px contours. Will need to be 3px bigger in w x h  e.g. 22x22 for 19x19 strel_sz ROIs
        SE = strel('square', strel_sz(1)+3); %<----------------------
        Afr=imclose(Afr,SE);  %<----------------------
        
        %     fr = 1111; %for TESTING
        %--Label connected components--
        BW = im2bw(Afr);  %make binary
        [L, num] = bwlabel(BW, 8);  %label connected regions
        RGB = label2rgb(L); %for display testing
        %     figure; imshow(RGB); %for display TESTING
        
        %--Get connected components areas--
        STATS = regionprops(L, 'Area'); %find areas
        areas=[STATS.Area]; %areas as single vector
        
        %--Delete connected components that are single pixels (1 strel_sz with zero neighbors)
        tmp = Afr;
        minIdx = find(areas <= 1*prod(strel_sz+1));  %get labeled region indices that are less than or equal to one strel_sz area with the 1 pixel border added by imdilate above
        if ~isempty(minIdx)
            for minIdx1 = minIdx
                tmp(L==minIdx1) = 0;
            end
        end
        %     figure; imshow(tmp);
        Afr = tmp;
        
        %----dilate connected regions to help with centroid calculations---
        %             SE = strel('disk',10);
        SE = strel('disk',round(strel_sz(1)*2.5));  %2.5x bigger than strel_sz. i.e. for 4px square strel that give approx 20um diam regions, this will give 10px diam that will give a 50-80um region.
        A3fr=imdilate(Afr,SE);
        SE = strel('disk',2);
        A3fr=imclose(A3fr,SE);
        %             figure;imshow(A2(:,:,fr))  %TESTING
        
        A2(:,:,fr) = Afr;
        A3(:,:,fr) = A3fr;
        
    end
    
    if rem(fr,10) == 0
        waitbar(fr/szZ,hbar);
    end
end
close(hbar)