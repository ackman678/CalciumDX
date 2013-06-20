%JBA, Tuesday, July 21, 2009 6:17 PM
%file_list = clipboard('pastespecial')  %if you have copied filelist to clipboard, array of an array
%file_list = file_list.A_pastespecial;

%file_list = clipboard('pastespecial'); file_list = file_list.A_pastespecial;
%Or if you have a file list as a plain text file you can use the import data
%graphic wizard

%once you have your file_list cell array created then perform the following
p_val = 0.05;
if exist('file_list','var') == 1
    if iscell(file_list)  %if multiple files are selected
        for i=1:length(file_list)
            fnm=file_list{i};
            load(fnm);
            
            if isfield(region,'transients') == 1
                region.userdata.active=find(region.transients > 1)';
                region.userdata.schmutzr=find(region.transients == 2)';
                region.userdata.gdp=find(region.transients == 3)';
                region.userdata.other=find(region.transients == 4)';
                region.userdata.sink=find(region.transients == 5)';
                region.userdata.schmutzon=cell(1,length(region.userdata.schmutzr));
                region.userdata.schmutzoff=cell(1,length(region.userdata.schmutzr));
                for c= 1:length(region.userdata.schmutzr)
                    region.userdata.schmutzon{c}=region.onsets{region.userdata.schmutzr(c)};
                    region.userdata.schmutzoff{c}=region.offsets{region.userdata.schmutzr(c)};
                end
            end
            numreshuffle = size(region.traces,2);
            corr_pairs = find_calciumdx_corr_pairs(region,1,numreshuffle,p_val)
            region.userdata.corr_pairs=corr_pairs;
            region.userdata.corrpval = p_val;
            
            save(fnm,'region')
        end
    end
    
else

    if exist('calciumdxprefs.mat','file') == 2
        load('calciumdxprefs')
    else
        pathname = pwd;
    end
    
    [filename pathname] = uigetfile('*.mat','Select calciumdx .mat file',pathname);
    fnm = fullfile(pathname,filename);
    load(fnm);
    save('calciumdxprefs.mat', 'pathname','filename');
    
    if isfield(region,'transients') == 1
        region.userdata.active=find(region.transients > 1)';
        region.userdata.schmutzr=find(region.transients == 2)';
        region.userdata.gdp=find(region.transients == 3)';
        region.userdata.other=find(region.transients == 4)';
        region.userdata.sink=find(region.transients == 5)';
        region.userdata.schmutzon=cell(1,length(region.userdata.schmutzr));
        region.userdata.schmutzoff=cell(1,length(region.userdata.schmutzr));
        for c= 1:length(region.userdata.schmutzr)
            region.userdata.schmutzon{c}=region.onsets{region.userdata.schmutzr(c)};
            region.userdata.schmutzoff{c}=region.offsets{region.userdata.schmutzr(c)};
        end
    end
    numreshuffle = size(region.traces,2);
    corr_pairs = find_calciumdx_corr_pairs(region,1,numreshuffle,p_val)
    region.userdata.corr_pairs=corr_pairs;
    region.userdata.corrpval = p_val;
    
    save(fnm,'region')
    
end


