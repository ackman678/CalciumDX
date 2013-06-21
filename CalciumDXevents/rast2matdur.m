function rast2matdur = rast2matdur(spk,endpt,sz)
%rast2mat = rast2mat(spk,endpt,dur)
%   converts a rasterplot into a matrix

s = zeros(length(spk),sz);
for c = 1:length(spk)
   for d = 1:length(spk{c})
      s(c,spk{c}(d):endpt{c}(d)) = 1;
   end
end

rast2matdur = s;