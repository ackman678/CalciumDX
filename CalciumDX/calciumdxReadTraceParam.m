    delete(bord_title)
    delete(regax)
    delete(det_tx1)
    delete(dpfilters)
    delete(det_loc)
    delete(det_view)
    delete(dummyp)
    delete(txlab)
    delete(txthres)
    delete(txarlow)
    delete(txarhigh)
    delete(btdetect)
    delete(bthide)
    delete(txpilim)
    delete(btfindbad)
    delete(btadjust)
    delete(btprev)
    delete(btnext)
    % delete(shaperad)
    delete(btadd)
    delete(btdelete)
    delete(btnextscr)
    delete(btclear)
    delete(btimport)
    delete(btalign)
    delete(btalign2)
    delete(btalign3)
    delete(btoggle1)
%    delete(btcalciumdxsavebackup)
%    delete(btcalciumdxloadbackup)

button = questdlg({'Divide ROIs into rectangular grid ROIs?'},'Rectangular ROI mesh','Yes','No','Yes');
if strcmp(button,'Yes')
    %use rectangular mesh
    region.cutoff = thres;
    region.lowarea = lowar;
    region.higharea = highar;
    region.isdetected = isdetected;
    region.pilimit = pilim;
    region.isadjusted = isadjust;
    region.contours = {};
    region.location = [];
    for c = 1:numel(cn)
        if isempty(cn{c})
            button2 = questdlg({['Make rectangular mesh grid over all of ' region.name{c} '?']},'Rectangular ROI mesh','Yes','No','Yes');
            if strcmp(button2,'Yes')
                cn{c}{1} = region.coords{c};
            end
        end
        if ~isempty(cn{c})
            for d= 1:numel(cn{c})
                %             if polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) > lowar(c) & polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) < highar(c)
                region.coords{c} = cn{c}{d};
                
                prompt = {'height','width'};
                dlg_title = 'size of rectangular ROI- default 20x40px for 256x512px,2xzoom; 10x20px for 256x512 1xzoom'; %10x20px for 1x zoom
                num_lines = 1;
                
                [szX,szY] = size(region.image);  %assuming szY is the largest dimension
                if mod(max([szY szX]),min([szY szX])) == 0
                    rXY=szY/szX;
                else
                    rXY = 1;
                end
                
                def_1=round(22/region.spaceres);
                def_1=def_1/rXY;
                def_2=round(22/region.spaceres);
                
                %         def = {'20','40'};
                def = {num2str(def_1),num2str(def_2)};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                strel_sz= [str2double(answer{1}) str2double(answer{2})];
                
                
                ps = round(cn{c}{d});
                [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
                inp = inpolygon(subx,suby,cn{c}{d}(:,1),cn{c}{d}(:,2));
                fx = subx(inp==1);
                fy = suby(inp==1);
                f = sub2ind(size(region.image),fy,fx);
                
                B=zeros(size(region.image));
                B(f)=1;
%				figure; imshow(B) %testing
                
                B2 = bwperim(B);
                [I,J]=ind2sub(sz,find(B2==1));
%				figure; imshow(B2) %testing
%		 hold on %testing
                
                I2=min(I):strel_sz(1):max(I);   %here is the problem
                J2=min(J):strel_sz(2):max(J);   %here is the problem	
                
                I2arr=repmat(I2',1,numel(J2));
                J2arr=repmat(J2,numel(I2),1);
                
                Arr=zeros(sz);
                Arr(I2arr,J2arr)=1;
                
                A = Arr .* B;   %much simpler for finding the interesecting set!
%				figure; imshow(A);                
                
                %the following is the old algorithm-- worked-- but only for convex polygons, not concave polygons (like crescent shapes)
%{
                [i,j]=ind2sub(sz,find(Arr==1));
                
                K=inpolygon(i,j,I,J);  %find those that are bounded by the original ROI image mask
         plot(j(K),i(K),'.','color','white') %shows the pixels bounded by the mask.  %testing
                
                A=zeros(sz);
                idx=sub2ind(sz,i(K),j(K));
                A(idx)=1;
                
%                         figure; imshow(A)  %testing
                %         % A(B2)=1;
                %         % figure; imshow(A)
%}
                
                % SE = strel('octagon', 9);
                SE = strel('rectangle', [strel_sz(1)-1 strel_sz(2)-1]);
                
                A2=imdilate(A,SE);
                figure;imshow(A2)
                
                [L,num]=bwlabel(A2);
                RGB = label2rgb(L);
                %         figure; imshow(RGB)
                
                
                %                 figure; subplot(2,1,1)
                %                 imshow(mat2gray(region.image))
                %                 hold on
                %                 plot(J,I,'.','color','red')
                %                 subplot(2,1,2)
                %                 imshow(RGB)
                %                 hold on
                %                 plot(J,I,'.','color','red')
                
                %get subindices row,col vectors for ROI no. 2
                BWall=zeros(sz);
                for i=1:num
                    %             region.location = [region.location d]; %numeric index for which region the contour is in.  Usually a vector of 1's the length of region contours
                    region.location = [region.location c];
                    [rows, cols] = find(L==i);
                    %             region.contours{length(region.contours)+1} = [rows cols];
                    
                    Atmp=zeros(sz);
                    Atmp(rows,cols)=1;
                    %             Atmp(region.contours{i}(:,1),region.contours{i}(:,2))=1;
                    %             Atmp(region.contours{i}(:,2),region.contours{i}(:,1))=1;
                    %             BWtmp=bwperim(Atmp);
                    %             [rows,cols]=ind2sub(sz,find(BWtmp));
                    %             region.contours{length(region.contours)+1} = [cols rows];
                    
                    %             figure;imshow(Atmp)
                    BWtmp = bwboundaries(Atmp,'noholes');
                    region.contours{length(region.contours)+1} = [BWtmp{1}(:,2) BWtmp{1}(:,1)];
                    %             hold on
                    %             plot(BWtmp{1}(:,2), BWtmp{1}(:,1), 'r', 'LineWidth', 2)
                    
                    %             BWall(BWtmp)=1;
                end
                %         figure; imshow(BWall)
                
                
                %             end
                
            end
        else
            continue
        end
    end
    
    halo_hands = [];
    halos = [];
    region.halomode = 0;
    region.haloarea = 1;
    
    trace_title = uicontrol(fig,'Style','text','Units','normalized','String','Traces','Position',[.87 .755 .11 0.03],'FontSize',12, ...
        'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
    
    halo_check = uicontrol(fig,'Style','checkbox','Units','normalized','String','Use halos','Position',[.87 .715 .11 0.025],'FontSize',9,...
        'BackgroundColor',[.8 .8 .8],'Callback','calciumdxHaloCheck');
    dummy(1) = uicontrol(fig,'Style','text','Units','normalized','String','Halo area','Position',[.87 0.6875 .11 0.02],'FontSize',9,...
        'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
    inpthaloar = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(region.haloarea),'Position',[.87 0.6625 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left','enable','off');
    btupdate = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Update','Position',[.93 .6175 .05 0.03],'FontSize',9,...
        'Callback','calciumdxHaloUpdate','enable','off');
    
    currdir = pwd;
    cd(fullfile(calciumdxpath,'TraceReaders'));
    mt = dir('*.m');
    cd(currdir);
    
    st = cell(1,length(mt));
    for c = 1:length(mt)
        st{c} = mt(c).name(1:end-2);
        if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'calciumdxTR_')
            st{c} = st{c}(9:end);
        end
    end
    
    dummy(2) = uicontrol(fig,'Style','text','Units','normalized','String','Trace reader','Position',[.87 0.575 .11 0.02],'FontSize',9,...
        'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
    dpreaders = uicontrol(fig,'Style','popupmenu','Units','normalized','String',st,'Position',[.87 .55 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1]);
    
    bnext = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 .03],'FontSize',9,'Callback','calciumdxReadTraces');
    
elseif strcmp(button,'No')
    %don't use rectangular mesh, use the manually defined ROIs
    
    if isimported < 1
        region.cutoff = thres;
        region.lowarea = lowar;
        region.higharea = highar;
        region.isdetected = isdetected;
        region.pilimit = pilim;
        region.isadjusted = isadjust;
        region.contours = {};
        region.location = [];
        for c = 1:length(cn)
            for d = 1:length(cn{c})
                if polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) > lowar(c) && polyarea(cn{c}{d}(:,1),cn{c}{d}(:,2)) < highar(c)
                    region.contours{length(region.contours)+1} = cn{c}{d};
                    region.location = [region.location c];
                end
            end
        end
    else
        disp('Using the region.coords, contours, and location info imported from file')
        region.contours = {};
        region.location = [];
        for c = 1:length(cn)
            for d = 1:length(cn{c})
                    region.contours{length(region.contours)+1} = cn{c}{d};
                    region.location = [region.location c];
            end
        end
        
    end
    
    halo_hands = [];
    halos = [];
    region.halomode = 0;
    region.haloarea = 1;
    
    trace_title = uicontrol(fig,'Style','text','Units','normalized','String','Traces','Position',[.87 .755 .11 0.03],'FontSize',12, ...
        'FontWeight','Bold','BackgroundColor',[.8 .8 .8]);
    
    halo_check = uicontrol(fig,'Style','checkbox','Units','normalized','String','Use halos','Position',[.87 .715 .11 0.025],'FontSize',9,...
        'BackgroundColor',[.8 .8 .8],'Callback','calciumdxHaloCheck');
    dummy(1) = uicontrol(fig,'Style','text','Units','normalized','String','Halo area','Position',[.87 0.6875 .11 0.02],'FontSize',9,...
        'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
    inpthaloar = uicontrol(fig,'Style','edit','Units','normalized','String',num2str(region.haloarea),'Position',[.87 0.6625 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left','enable','off');
    btupdate = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Update','Position',[.93 .6175 .05 0.03],'FontSize',9,...
        'Callback','calciumdxHaloUpdate','enable','off');
    
    currdir = pwd;
    cd(fullfile(calciumdxpath,'TraceReaders'));
    mt = dir('*.m');
    cd(currdir);
    
    st = cell(1,length(mt));
    for c = 1:length(mt)
        st{c} = mt(c).name(1:end-2);
        if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'calciumdxTR_')
            st{c} = st{c}(9:end);
        end
    end
    
    dummy(2) = uicontrol(fig,'Style','text','Units','normalized','String','Trace reader','Position',[.87 0.575 .11 0.02],'FontSize',9,...
        'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
    dpreaders = uicontrol(fig,'Style','popupmenu','Units','normalized','String',st,'Position',[.87 .55 .11 0.025],'FontSize',9,...
        'BackgroundColor',[1 1 1]);
    
    bnext = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next >>','Position',[.93 .02 .05 .03],'FontSize',9,'Callback','calciumdxReadTraces');
end