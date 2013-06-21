function rast2mat = rast2mat(spk,sz)
%rast2mat = rast2mat(spk)
%   converts a rasterplot into a matrix

s = zeros(size(spk,2),sz);
for c = 1:size(spk,2)
   s(c,spk{c}) = 1;
end

rast2mat = s;