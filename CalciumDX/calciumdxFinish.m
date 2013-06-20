[filename, pathname] = uiputfile('*.mat', 'Save file as',fnm(1:end-4));
if ~ischar(filename)
    return
end
fnm = [pathname filename];

save temp.mat region fnm pathname filename
set(gcf,'CloseRequestFcn','');
delete(gcf)
clear
load temp.mat region fnm pathname filename
delete temp.mat
save(fnm,'region');
clear fnm

load('calciumdxbackup.mat')
delete('calciumdxbackup.mat')


load calciumdxprefs
if isdir('../CalciumDXevents')
    cd ../CalciumDXevents
    calciumdxevents
elseif isdir('./CalciumDXevents')
    cd ../CalciumDXevents
    calciumdxevents
else
    error('calciumdxevents folder not found')
end

