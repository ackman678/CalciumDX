function [s_reorder] = rasterReorder(s,revalueSpk)
%---convert bin raster from vertical to horizontal roi labeling to see if structural patterns in raster has improved clarity
%James B. Ackman 2011-12-28
if nargin < 2 || isempty(revalueSpk), revalueSpk = 0; end
s_reorder = zeros(size(s));
k = 0;
spikeValue=0;
for j = 1:10
    ind = j:10:100;
    for i = 1:numel(ind)
        k = k+1;
        s_reorder(k,:) = s(ind(i),:);
        
        if revalueSpk > 0
        spikeValue = spikeValue+1;
        [m n] = find(s_reorder(k,:));
        s_reorder(k,n) = spikeValue;
        end
    end
end

k = 100;
spikeValue=0;
for j = 101:110
    ind = j:10:200;
    for i = 1:numel(ind)
        k = k+1;
        s_reorder(k,:) = s(ind(i),:);
        
        if revalueSpk > 0
        spikeValue = spikeValue+1;
        [m n] = find(s_reorder(k,:));
        s_reorder(k,n) = spikeValue;
        end
    end
end

I1 = s_reorder(1:100,:);
I1 = flipud(I1);
I2 = s_reorder(101:200,:);
% I2 = flipud(I2);
s_reorder = [I1; I2];