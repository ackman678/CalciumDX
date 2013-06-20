figure; imagesc(region.traces)

figure; 
a(1)=subplot(2,1,1);
plot(nt(68,:))
a(2)=subplot(2,1,2);
plot([diff(nt(68,:)) 0])
linkaxes(a,'x')


A = [];
for i=1:size(region.traces,1)
for c=1:length(region.onsets{i})
X = region.onsets{i}(c):region.offsets{i}(c);
Y = nt(i,X);
X = [X X(end) X(1)];
Y = [Y Y(1) Y(1)];
% figure; plot(X,Y)
% A = polyarea(X,Y);
A = [A polyarea(X,Y)];
end
end
figure; hist(A,50)

tmp = trHist;
waveframes = [];
if ~isempty(artifactframes) && isempty(waveframes)    
    waveframes = setxor(artifactframes(:,1), [1:length(region.traces)])';
    c=intersect(find(trHist > 0),waveframes(:,1)); %this will limit the selected possible waveframes, to just those that could have occured during a wave (those frames with 1 or more active ROIs).
    d = setxor(c, [1:length(region.traces)]);
if ishandle(hf1(1))
    close(hf1(1))
end
    trHist(d)=0;
    figure;
    a(1) = subplot(2,1,1)
    plot(trHist)
    a(2) = subplot(2,1,2)
    imagesc(nt)
    linkaxes(a,'x')
end


if isempty(artifactframes) && ~isempty(waveframes)
    c=intersect(find(trHist > 0),waveframes(:,1)); %this will limit the selected possible waveframes, to just those that could have occured during a wave (those frames with 1 or more active ROIs).
    d = setxor(c, [1:length(region.traces)]);
if ishandle(hf1(1))
    close(hf1(1))
end
    trHist(d)=0;
    figure;
    a(1) = subplot(2,1,1)
    plot(trHist)
    a(2) = subplot(2,1,2)
    imagesc(nt)
    linkaxes(a,'x')
end




%-----------Try image filter on imagesc of traces for artifact detection-----------------------------
A = nt;
SE = strel('rectangle', [fix((size(nt,1))) 1]);
A2 = imclose(A,SE);
% SE = strel('ball', 10, 3);
% A2=imdilate(A,SE);
% figure;imshow(A2)
% A2=imdilate(A,SE);
% SE = strel('rectangle', [5 size(nt,2)]);
% A2=imerode(A2',SE);
figure;
% a(1) = subplot(3,1,1)
% imagesc(nt)
% a(2) = subplot(3,1,2)
a(1) = subplot(2,1,1)
imagesc(A2)
% a(3) = subplot(3,1,3)
a(2) = subplot(2,1,2)
I = mat2gray(A2);
BW = edge(I);
imshow(BW)

F = fft2(A);
Fc = fftshift(F);
S2 = log(1 + abs(Fc));
% figure; imagesc(abs(Fc))
figure; imagesc(S2)
% figure; imshow(abs(Fc),[])

%------fft notch filter code-----------------------------------------------
A = nt;
% A = region.traces;
% [M,N] = size(f);
[M,N] = size(A);
% [f,revertclass] = tofloat(f);
% [f,revertclass] = tofloat(A);
% F = fft2(f);
F = fft2(A);
S = gscale(log(1+abs(fftshift(F))));
% figure; imshow(S);  %show fft image
% H = recnotch('reject','horizontal',M,N,21,15,15);
H = recnotch('reject','horizontal',M,N,21,31,31);
% figure; imshow(fftshift(H))  %show what the notch filter on the fft image will be
%display the filtered image
g = dftfilt(A,H);
% H = recnotch('reject','vertical',M,N,11,31,31);
% g = dftfilt(g,H);
% g = revertclass(g);
% figure; imshow(g);
figure;
a(1) = subplot(2,1,1);
imagesc(A)
a(2) = subplot(2,1,2);
imagesc(g);
linkaxes(a,'x')
nt = g;
region.tracesRaw = region.traces;
region.tracesFilt = nt;
%--------------------------------------------------------------------------

%pass the interference image instead
Hrecpass = recnotch('pass', 'horizontal', M, N, 21, 31,31);
interference = dftfilt(A, Hrecpass);
figure, imshow(fftshift(Hrecpass))
interference = gscale(interference);
figure, imshow(interference);
%---END template code-----------------------------------------


%---template code-----------------------------------------
%use custom recnotch filter from Gonzalez et al. Digital Image Processing Book
[M,N] = size(f);
[f,revertclass] = tofloat(f);
F = fft2(f);
S = gscale(log(1+abs(fftshift(F))));
imshow(S);
H = recnotch('reject','vertical',M,N,3,15,15);
figure; imshow(fftshift(H))
%display the filtered image
g = dftfilt(f,H);
g = revertclass(g);
figure; imshow(g);

%pass the interference image instead
Hrecpass = recnotch('pass', 'vertical', M, N, 3, 15,15);
interference = dftfilt(f, Hrecpass);
figure, imshow(fftshift(Hrecpass))
interference = gscale(interference);
figure, imshow(interference);
%---END template code-----------------------------------------



%---test wavelet decomposition---------------------------------
s = data;
% Decompose the signal s at level 5 using the wavelet db3. 
w = 'mexh'; 
[c,l] = wavedec(s,5,w);

% Reconstruct the details using the decomposition structure. 
for i = 1:5
   D(i,:) = wrcoef('d',c,l,w,i);
end

figure;
%Avoid edge effects by suppressing edge values and plot. 
tt = 1+100:length(s)-100; 
subplot(6,1,1); plot(tt,s(tt),'r'); 
title('Electrical Signal and Details'); 
for i = 1:5, subplot(6,1,i+1); plot(tt,D(5-i+1,tt),'g'); end
%-----------------------------------