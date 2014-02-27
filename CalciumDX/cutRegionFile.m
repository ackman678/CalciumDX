function region = cutRegionFile(frames,filenameIN,filenameOUT)
%cutRegionFile - Batch processing script for subsetting a region data structure
%Examples:
% >> cutRegionFile;
% >> cutRegionFile([1 150]);
% >> region = cutRegionFile([1 150]);
% >> cutRegionFile([350 700],'file1.mat','file2_fr350-700.mat')
%
%**USE**
% Options: 
% frames - a two element vector containing the first frame and last frame you want to split your movie into.
% filenameIN - string, full filename of the .mat file containing the region data structure you want to cut
% filenameOUT - string, full filename of the .mat file you want to save the cut region data structure as
%James B. Ackman 2014-02-27 13:50:46

if nargin < 1 || isempty(frames),
	prompt = {'first frame:','last frame:'};
	dlg_title = 'Input new frame indices';
	num_lines = 1;

	def = {'1','150'};

	optionsDlg.Resize='on';
	answer = inputdlg(prompt,dlg_title,num_lines,def,optionsDlg);

	frames(1) = str2num(answer{1});
	frames(2) = str2num(answer{2});
end
%{
[filename, pathname] = uigetfile('*.mat','MultiSelect','on');
if iscell(filename)  %if multiple files are selected
	multiFile(filename)	
	for i=1:length(filename)
		fnm = [pathname filename{i}];
		processFile(fnm,frames);
	end	
elseif ischar(filename)   %if only one file is selected
	fnm = [pathname filename];
	processFile(fnm,frames);
else 
	error('wrong input')
end
%}

if nargin < 2 || isempty(filenameIN)
	[filename, pathname] = uigetfile('*.mat');
	fnm = fullfile(pathname,filename);
else
	fnm = filenameIN;
end

if nargin < 3 || isempty(filenameOUT)
	fnm2 = [fnm(1:end-4) '_fr' num2str(frames(1)) '-' num2str(frames(2)) '.mat'];
else
	fnm2 = filenameOUT;
end

if ischar(filename)   %if only one file is selected
	region = processFile(fnm,frames);
else 
	error('wrong input')
end

save(fnm2,'region')


function region = processFile(fnm,frames)
load(fnm);
sz = size(region.traces);
if frames(1) < 1 || frames(1) > sz(2) || frames(2) < 1 || frames(2) > sz(2)
	error('New frame indices not within range (must be between 1 and size(region.traces,2))')
end

region.traces = region.traces(:,frames(1):frames(2));

newIndices = frames(1):frames(2);
indShift = frames(1) - 1; 

if ~isempty(region.onsets)  %determine whether region.onsets contains data (whether signal detection has been done yet)
	for c = 1:size(region.onsets,2) %fix onset and offset values
		disp('============================================')
		disp(c)
		[idxMatch1,ia1,ib1] = intersect(newIndices,region.onsets{c});
		[idxMatch2,ia2,ib2] = intersect(newIndices,region.offsets{c});

		if ~isempty(idxMatch1) && (isempty(idxMatch2) | length(idxMatch2) < length(idxMatch1))  %A single on-off overlapping at end of section
			region.offsets{c} = [region.offsets{c} frames(2)];
			[idxMatch1,ia1,ib1] = intersect(newIndices,region.onsets{c});
			[idxMatch2,ia2,ib2] = intersect(newIndices,region.offsets{c});
			disp('End overlap fix')
		end
		
		if (isempty(idxMatch1) | length(idxMatch1) < length(idxMatch2))  && ~isempty(idxMatch2)  %A single on-off overlapping at start of section
			region.onsets{c} = [frames(1) region.onsets{c}];
			[idxMatch1,ia1,ib1] = intersect(newIndices,region.onsets{c});
			[idxMatch2,ia2,ib2] = intersect(newIndices,region.offsets{c});
			disp('Start overlap fix')
		end

		if ~isempty(idxMatch1) && ~isempty(idxMatch2)
%			if ib2(1) < ib1(1)  %if a transient on-off is overlapping at start of section
%				region.onsets{c} = [frames(1) region.onsets{c}(ib1)];  %TODO: Note this might set the onset and offset to have same value
%				region.offset{c} = region.offsets{c}(ib2);
%				disp('third overlap fix')
%			else
%				region.onsets{c} = region.onsets{c}(ib1);
%				region.offsets{c} = region.offsets{c}(ib2);
%			end
%			if ib2(end) < ib1(end) %if a transient on-off is overlapping at end of section
%				region.onsets{c} = region.onsets{c}(ib1);
%				region.offsets{c} = [region.offsets{c}(ib2) frames(2)];  %TODO: Note this might set the onset and offset to have same value
%				disp('fourth overlap fix')
%			else
%				region.onsets{c} = region.onsets{c}(ib1);
%				region.offsets{c} = region.offsets{c}(ib2);		
%			end
%			region.onsets{c} = region.onsets{c} - indShift;
%			region.offsets{c} = region.offsets{c} - indShift;

			region.onsets{c} = region.onsets{c}(ib1);
			region.offsets{c} = region.offsets{c}(ib2);
		
			region.onsets{c} = region.onsets{c} - indShift;
			region.offsets{c} = region.offsets{c} - indShift;

		else
			
			region.onsets{c} = [];
			region.offsets{c} = [];
		end
	end
end
