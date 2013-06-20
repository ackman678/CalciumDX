set(bord_delete,'foreground',[1 0 0]);
[x y butt] = ginput(1);
if butt > 1
    set(bord_delete,'foreground',[0 0 0]);
    return
end

dst = [];
for c = 1:length(bord)
    dst(c) = min(sum((bord{c}-repmat([x y],size(bord{c},1),1)).^2,2));
end
[mn i] = min(dst);

delete(bhand(i));
bhand = bhand([1:i-1 i+1:length(bord)]);

tmp = [];
for c = [1:i-1 i+1:length(bord)]
    tmp{length(tmp)+1} = bord{c};
end
bord = [];
for c = 1:length(tmp)
    bord{c} = tmp{c};
end

set(bord_delete,'foreground',[0 0 0]);
calciumdxDetermineRegions