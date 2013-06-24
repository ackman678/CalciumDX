[filename, pathname] = uiputfile('*.mat', 'Save file as',fnm(1:end-4));
if ~ischar(filename)
    return
end
fnm = [pathname filename];
save(fnm,'region');
save(calciumdxprefs,'pathname', 'filename')

save(calciumdxbackupLocation, 'region', 'fnm', 'pathname', 'filename')
set(gcf,'CloseRequestFcn','');
delete(gcf)
clear

matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxbackupLocation = fullfile(matlabUserPath,'calciumdxbackup.mat');

load(calciumdxbackupLocation)
delete(calciumdxbackupLocation)

%if isdir(fullfile('..','CalciumDXevents'))
%    cd(fullfile('..','CalciumDXevents'))
%    calciumdxevents
%elseif isdir(fullfile('.','CalciumDXevents'))
%    cd(fullfile('.','CalciumDXevents'))
%    calciumdxevents
%else
%    error('calciumdxevents folder not found')
%end

cd(matlabUserPath)
clear calciumdxbackupLocation matlabUserPath
calciumdxevents
