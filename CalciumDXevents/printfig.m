function printfig(type,fname,homedir,hires)
%printfig
%Figure export and save wrapped around the print function.
%type -- can be any of '-depsc', '-djpg', '-dpng', '-dtiff', '-dtiffn', '-dsvg, '-dpdf', or any others as described in help(print)
%fname -- optional provided filename
%homedir -- 'true' or 'false', whether or not to save to users' home matlab directory (userpath variable). Defaults to 'true'
%hires -- whether to use high resolution for png images. Defaults to 'false'
%
%James B. Ackman Monday, November 28, 2011 5:43 PM
%Fixed homedir option for all OSs and added hires as an option 2013-07-17 17:53:56

if nargin < 4 || isempty(hires), hires = 'false'; end
if nargin < 3 || isempty(homedir), homedir = 'true'; end
if nargin < 2 || isempty(fname), fname = [datestr(now,'yyyymmdd-HHMMSS') '_figure' num2str(gcf)]; end
if nargin < 1 || isempty(type), type = 'epsc'; end
if strcmp(homedir,'true')
    [pathstr,name,ext]=fileparts(fname);
    fname = [name ext];
    matlabUserPath = userpath;
    matlabUserPath = matlabUserPath(1:end-1);
    fname = fullfile(matlabUserPath,fname);
end
% saveas(gcf,fname, type);
% print('-depsc','-tiff','-r300',fname)
% print(gcf, '-dpng', '-r0', fname);
set(gcf,'PaperPositionMode','auto')
if strcmp(type,'epsc')
    print(gcf, '-depsc', fname);
elseif strcmp(type,'png') && strcmp(hires,'true')
    print(gcf, '-dpng', '-r600', fname);
elseif strcmp(type,'png')
    print(gcf, '-dpng', fname);
else
    print(gcf, type, '-r0', fname);
end
