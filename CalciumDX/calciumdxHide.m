if ishid == 0
    set(bthide,'string','Show');
    set(handl{num},'visible','off');
else
    set(bthide,'string','Hide');
    calciumdxDrawCells
end