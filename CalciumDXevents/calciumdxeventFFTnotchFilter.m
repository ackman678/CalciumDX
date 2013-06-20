%calciumdxeventFFTnotchFilter.m
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
figure; imshow(S);  %show fft image
% H = recnotch('reject','horizontal',M,N,21,15,15);
H = recnotch('reject','horizontal',M,N,21,31,31);
figure; imshow(fftshift(H))  %show what the notch filter on the fft image will be
%display the filtered image
g = dftfilt(A,H);
% H = recnotch('reject','vertical',M,N,11,31,31);
% g = dftfilt(g,H);
% g = revertclass(g);
figure; imshow(g);
figure;
a(1) = subplot(2,1,1);
imagesc(A)
a(2) = subplot(2,1,2);
imagesc(g);
linkaxes(a,'x')
nt = double(g);
region.tracesRaw = region.traces;
region.tracesFilt = nt;
%--------------------------------------------------------------------------