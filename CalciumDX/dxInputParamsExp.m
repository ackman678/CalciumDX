prompt = {'Type (wt, C57BL6, GAD67.GFP, etc):','Age (P7, adult, etc):','Type of experiment (2P.img, CCD.img, or img.ephys):','Calcium dye used:','Brain area (SC.L, V1.R, ctx, etc):','Enter field factor (f1, f2, etc):','Enter y or n if z-artifact present:','Enter z-depth below pia (in um):','Enter type of anesthetic:','Enter inhal anesthetic percentage:','Comments:'};
dlg_title = 'Input experimental parameters';
num_lines = 1;
def = {'wt','','img','CaGreenDextran','SC.R','','n','','isoflurane','1.0',''};
optionsDlg.Resize='on';
answer = inputdlg(prompt,dlg_title,num_lines,def,optionsDlg);

region.animaltype = answer{1};
region.age= answer{2};
region.exptype= answer{3};
region.dye= answer{4};
region.brainarea = answer{5};
region.field= answer{6};
region.zartifact= answer{7};
region.zdepth= answer{8};
region.anesthetic= answer{9};
region.anesthpercent= str2double(answer{10});
region.comments= answer{11};