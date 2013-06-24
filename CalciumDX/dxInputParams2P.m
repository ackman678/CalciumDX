prompt = {'Enter whether frames were averaged:','Enter framerate (Hz):','Enter no. scanned lines per frame:','Enter scan line period:','Enter optical zoom:'};
dlg_title = 'Input 2P experimental parameters';
num_lines = 1;
def = {region.frameaveraging,num2str(region.framerate),num2str(region.linesperframe),num2str(region.scanlineperiod),num2str(region.opticalzoom)};
optionsDlg.Resize='on';
answer = inputdlg(prompt,dlg_title,num_lines,def,optionsDlg);

region.frameaveraging= answer{1};
region.framerate= str2double(answer{2});
region.linesperframe= str2double(answer{3});
region.scanlineperiod= str2double(answer{4});
region.opticalzoom= str2double(answer{5});