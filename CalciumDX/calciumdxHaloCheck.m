if get(halo_check,'value') == 0
    delete(halo_hands);
    halo_hands = [];
    halos = [];
    region.halomode = 0;
    set(inpthaloar,'enable','off');
    set(btupdate,'enable','off');
else
    region.halomode = 1;
    set(btupdate,'enable','on');
    set(inpthaloar,'enable','on');
end
zoom on