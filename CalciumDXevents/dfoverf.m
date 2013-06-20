function dfoverf = dfoverf(tr)

bases = repmat(mean(tr,2),1,size(tr,2));
dfoverf = (tr-bases)./(bases+eps);

% if mean(f) == 0
%    dfoverf = f / (max(f)-min(f)+eps);
% else
%    a = (f - median(f))/median(f);
%    dfoverf = a ;
% end
