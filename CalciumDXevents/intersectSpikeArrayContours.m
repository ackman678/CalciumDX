function [s] = intersectSpikeArrayContours(contours,movie)
%James B. Ackman 2011-08-29
%to create new 2D spike array from 3D spike array movie, in effect spatially smoothing the resulting dataset
s = zeros(length(contours),size(movie,3));
% i = 25; fr=1111; movie = A2; contours = bigmesh.contours;    %TESTING
sz = size(movie);
% figure;  %TESTING
hbar = waitbar(0,'Please wait...');
% tic; inArr = find(movie); toc;
for i = 1:length(contours)
    ps = round(contours{i});
    [subx suby] = meshgrid(min(ps(:,1)):max(ps(:,1)),min(ps(:,2)):max(ps(:,2)));
    inp = inpolygon(subx,suby,ps(:,1),ps(:,2));
    fx = subx(inp==1);
    fy = suby(inp==1);
    f = sub2ind(sz(1:2),fy,fx);
    B=zeros(sz(1:2),'uint8');
    B(f)=1;
    
%     imshow(mat2gray(B))  %TESTING
%     drawnow  %TESTING

%     tic
    for fr = 1:sz(3)
        Afr = movie(:,:,fr);
        inB = find(B);
        inAfr = find(Afr);
        if ~isempty(intersect(inAfr,inB));
            s(i,fr) = 1;   %make into a counter instead for total no. of spikes within ea bigmesh.contour? Would this help xcorr value calculations? Would have to change the signal smoothing operations though.
        end
        clear inB inAfr
    end
%     toc 
    
%     tic; Barr = repmat(B,[1 1 sz(3)]); toc;
%     tic; inB = find(Barr); toc;
%          
%     c = intersect(inArr,inB);
%     if ~isempty(c)
%         [I,J,K] = ind2sub(size(movie),c);
%         s(i,K) = 1;
%     end
%     clear inB I J K c
%     if rem(i,10) == 0
        waitbar(i/size(s,1),hbar);
%     end
end
close(hbar)  
% figure; imagesc(s);


% figure; 
% for i=1:length(contours)
%     plot(contours{i}(:,1),contours{i}(:,2),'-k')
%     hold on
% end
% axis ij
    