matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxprefs = fullfile(matlabUserPath,'calciumdxprefs.mat');

if exist('pathname','var')
    fnm = [pathname filename];
    [filename, pathname] = uiputfile('*.mat', 'Save file as',fnm(1:end-4));
    if ~ischar(filename)
%         delete(gcf)
%         clear
        disp('no filename given...')
    end
    if exist('calciumdxprefs.mat','file') == 2, save(calciumdxprefs,'pathname','filename','-append'); else, save(calciumdxprefs,'pathname','filename'); end
else
    [filename, pathname] = uiputfile({'*.mat'}, 'Save file as');
    if ~ischar(filename)
%         delete(gcf)
%         clear
        disp('no filename given...')
    end
    if exist('calciumdxprefs.mat','file') == 2, save(calciumdxprefs,'pathname','filename','-append'); else, save(calciumdxprefs,'pathname','filename'); end
end

fnm = [pathname filename];


region.onsets = cell(1,size(spk,1));
region.offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    region.onsets{c} = find(spk(c,:)==1);
    region.offsets{c} = find(dec(c,:)==1);
end
save(fnm,'region');