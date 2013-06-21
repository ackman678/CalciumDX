NpopupDetect = get(popupDetect,'value');
detectorName = popupDetectList{NpopupDetect};
function_handle = str2func(detectorName);
hbar = waitbar(0,'Please wait...'); 
spk = zeros(size(nt));
dec = zeros(size(nt));
sz=size(region.traces);
for c = 1:sz(1) 
%     [s d] = calciumdxdettrialWaves(trSign*nt(c,:),region);
    [s,d,sd,sd2,sd3,nonfilt,hannfilterorder] = feval(function_handle,trSign*nt(c,:),region);
    if rem(c,10) == 0
        waitbar(c/sz(1),hbar); 
    end 
    spk(c,s) = 1; 
    dec(c,d) = 1; 
    if region.transients(1,c) == 1 && sum(spk(c,:)) > 0
        region.transients(1,c) = 4;
    end 
end
region.detectorname=['unsupervised calciumdxdettrial, sd=' num2str(sd) ', sd2=' num2str(sd2) ', sd3=' num2str(sd3) ', nonfilt=' num2str(nonfilt) ', HannOrder=' num2str(hannfilterorder)];
close(hbar)
hevPlotTrace
