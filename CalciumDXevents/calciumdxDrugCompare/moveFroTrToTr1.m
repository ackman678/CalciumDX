str = sum(get(gcf,'currentcharacter'));

if str == 28 | str == 44
    num1 = num1 - 1;
    set(txcellnum1,'string',num2str(num1));
    dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);
    xlimits = [0 size(nt1,2)+1];
    
    set(txcellnum2,'string',num2str(num1));
    dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);
    xlimits = [0 size(nt2,2)+1];
end

if str == 29 | str == 46
    num1 = num1 + 1;
    set(txcellnum1,'string',num2str(num1));
    dispCalcTracNetActComp(txcellnum1,region1,spk1,dec1,imgax1,trax1,nt1);
    xlimits = [0 size(nt1,2)+1];
    
    set(txcellnum2,'string',num2str(num1));
    dispCalcTracNetActComp(txcellnum2,region2,spk2,dec2,imgax2,trax2,nt2);
    xlimits = [0 size(nt2,2)+1];
end