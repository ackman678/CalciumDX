%JBA, Wednesday, October 22, 2008 2:41 PM
%returns a plot of the averaged Welch Power Spectral Density Estimate for a multiple channels (multiple calcium recorded cells)
%myCells is a vector of cell indices to select from your nt matrix of multple recordings (like active, non-active, gdp, spa, etc)
function PxxM = myPSD(nt,myCells,Fs)

if numel(myCells) > 0
    PxxM = [];
    for i=1:length(myCells)
%     [Pxx1, f] = pwelch(nt(myCells(i),:),hamming(200),[],1024,Fs);
    [Pxx1, f] = pwelch(nt(myCells(i),:),[],[],1024,Fs);
    %[Cxy1,F] = mscohere(r1,nt(myCells(i),:),[],[],1024,1/region.timeres);  %coherence of signal with gaussian noise
    PxxM(:,i) = Pxx1(:,1);   %add new values to matrix
%    h = waitbar(i/numel(myCells));
    end
    %assignin('base', 'PxxM',PxxM)
%    close(h)
    Hpsd = dspdata.psd(mean(PxxM,2),'Fs',Fs);   %mean of the Power spectrum values
    figure;
    plot(Hpsd)  %plot the mean power spectrum
else
    PxxM = [];
    print('Error-- less than one cell input')   
end
end