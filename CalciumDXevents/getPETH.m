function [sdf,Hist_raw,xout,meanLatencies,nSpikes,cells,frames,meanLatencyByCell]=getPETH(spk,AlignEv,bn,ktype,sig,binWidth,meanLatencies,nSpikes)
% function sdf=SDFConv(spk,AlignEv,bn,ktype,sig)
% 
% Inputs:
%   SPK: =  SPK is a 2D numeric array with the timestamp (in ms) for each 
%           spike in each trial down the other. SPK is assumed to be 
%           Nan-padded so that it's width is equal to the maximum number of
%           spikes in one trial.
%   ALIGNEV: =  A vector of the time of the Align event in each trial. The
%               program converts it to a column vector if it isn't already.
%               Defaults to all zeros (assuming the input is already
%               aligned)
%   BN:   = The time bin around the Align event you you would like
%           to analyze. In form [starttime endtime]. Defaults to [-300 800]
%           *Note that the effects from spikes from before bn(1) are 
%           included in the smoothed histogram, i.e. sdf, (including 
%           effects from spikes up to the kernel's half-width before 
%           bn(1)), though the time is not included in the output of the 
%           smoothed histogram
%   KTYPE: =    The type of kernel to use. 1 for PSP kernel, 2 for Gaussian
%               kernel. Defaults to 1
%   SIG: =  If using a gaussian kernel, this is the value of sigma (i.e.
%           the std of the kernel) in units of ms. Defaults to 30
%   
% OUTPUTS
%   SDF: =  The smoothed histogram
%
% Note- this doesn't require a for loop, and such easy translation into
% histograms is the main advantage to storing spikes in a 2D trials array
% padded with Nans, which is why it was determined to store spikes in that 
% way.
%
% Original SDFConv.m By Matt Nelson
%
% xout time bins added by James Ackman 7/22/2011 as well as hist raw and N's and means output
if nargin<8 || isempty(nSpikes);     nSpikes=[];              end
if nargin<7 || isempty(meanLatencies);     meanLatencies=[];              end
if nargin<6 || isempty(binWidth);     binWidth=10;              end
if nargin<5 || isempty(sig);     sig=30;              end
% if nargin<4 || isempty(ktype);   ktype=1;            end
if nargin<4 || isempty(ktype);   ktype=[];            end
if nargin<3 || isempty(bn);      bn=[-300 800];      end
if nargin<2 || isempty(AlignEv);      AlignEv=repmat(0,size(spk,1),1);      end

if ~isempty(ktype)
% Make the Kernel for the convolution.
if ktype==1     %psp kernel
    Growth=1; %ms
    Decay=20; %ms, to match grav PSP
    BW=round(Decay*8);
    
    Kernel=[0:BW];
    Kernel=(1-(exp(-(Kernel./Growth)))).*(exp(-(Kernel./Decay)));
    Kernel=Kernel./sum(Kernel);  %this makes the Kernel sum to 1
    Kernel=Kernel.*1000;    %this makes the Kernel sum to 1000... because the time bin of 1 msec has a time unit of 1 in our plots, this makes the output sdf in unit of Hz (spks/sec)
    
    edgedif=[BW 0];
else    %gaussian kernel
    BW=160;    %half the width of the gaussian kernel
    Kernel = normpdf([-BW:BW],0,sig);
    Kernel=Kernel./sum(Kernel).*1000;
    
    edgedif=[BW BW];
end
hbn=[bn(1)-edgedif(1) bn(2)+edgedif(2)];
else
hbn = bn;
end

%subtract the AlignEv for each trial from every timestamp
if size(AlignEv,2)>size(AlignEv,1)      AlignEv=AlignEv';       end     %make AlignEv a column vector
spk=spk-repmat(AlignEv,1,size(spk,2));

inds=spk>=hbn(1) & spk <= hbn(2);
[cells,frames]=find(spk>=hbn(1) & spk <= hbn(2));
[Hist_raw,xout]=hist(spk(inds),[hbn(1):binWidth:hbn(2)]);     %outputs 1 bin for ms
% figure; hist(spk(inds),[hbn(1):binWidth:hbn(2)])
sdf = Hist_raw;
meanLatencies = [meanLatencies; mean(spk(inds)) ];
nSpikes = [nSpikes; numel(spk(inds))];
meanLatencyByCell = min(spk(:,unique(frames)),[],2);
% disp(['mean latency = ' num2str(meanLatencies(end)) ' msec'])
% disp([ 'N spikes = ' num2str(nSpikes(end))])
% disp([ 'length hist = ' num2str(numel(Hist_raw))])

% sdf=filter(Kernel,1,Hist_raw) /size(spk,1);  % Convolve the raw histogram for each trial with the kernel + normalize by # of trials
% if ktype==1 
%     sdf=sdf(1+edgedif(1):end-edgedif(2));    %Cut off the 'breadcrusts'
% else
%     sdf=sdf(1+BW*2:end);   %Cut off the 'breadcrusts', and adjust timeshift for the Gaussian kernel
% end

