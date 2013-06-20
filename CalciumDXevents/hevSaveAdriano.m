[filename2, pathname2] = uiputfile([filename(1:end-3) 'mat'], 'Save file as');
if ~isstr(filename2)
    return
end
fnm = [pathname2 filename2];
filename = filename2;

region.contours=region.contours(setdiff(1:size(spk,1),cell2Del));
spk(cell2Del,:)=[];
dec(cell2Del,:)=[];
region.traces(cell2Del,:)=[];
region.location(cell2Del)=[];

region.onsets = cell(1,size(spk,1));
region.offsets = cell(1,size(dec,1));
for c = 1:size(spk,1)
    region.onsets{c} = find(spk(c,:)==1);
    region.offsets{c} = find(dec(c,:)==1);
end
save(fnm,'region');