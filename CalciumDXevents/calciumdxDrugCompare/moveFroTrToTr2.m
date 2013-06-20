str = sum(get(gcf,'currentcharacter'));

if str == 28 | str == 44
    num2 = num2 - 1;
    set(txcellnum1,'string',num2str(num2));
    dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);
    xlimits = [0 size(nt1,2)+1];
    
    set(txcellnum2,'string',num2str(num2));
    dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);
    xlimits = [0 size(nt2,2)+1];
end

if str == 29 | str == 46
    num2 = num2 + 1;
    set(txcellnum1,'string',num2str(num2));
    dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);
    xlimits = [0 size(nt1,2)+1];
    set(txcellnum2,'string',num2str(num2));
    dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);
    xlimits = [0 size(nt2,2)+1];
end