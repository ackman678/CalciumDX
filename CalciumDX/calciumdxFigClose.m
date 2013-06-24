clear
matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxbackupLocation = fullfile(matlabUserPath,'calciumdxbackup.mat');

load(calciumdxbackupLocation)
delete(gcf)
delete(calciumdxbackupLocation)

cd(matlabUserPath)
clear calciumdxbackupLocation matlabUserPath
