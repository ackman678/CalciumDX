function [sig_idx1 times] = myFrameTriggerDetect(fhandle,channelToFilter,framerate)
%Fetch frame times for a frame trigger pulse
%This frame trigger pulse should be an event channel in Spike2, for example giving times of frames from Prairie 2P or PCO CCD camera.
%[sig_idx1 times] = myFrameTriggerDetect(1,2,20);
%myPulseDetect.m
%James B. Ackman (c) 7/18/2011 

if nargin < 3, framerate=20; end %video frame rate you are trying to detect, to give a ballpark figure on the signal frequency we must detect

[fhandle, channels]=scParam(fhandle);

if ~isempty(channels{channelToFilter}.adc(:))
data=(channels{channelToFilter}.adc(:,1))';
% if nargin < 4; tolerance2=0.1; end
% if nargin < 3; tolerance1=0.1; end
stim_time = round(diff(data(1,:)));
% stim_time=stim_time*-1;

%data signal rate
Fs=getSampleRate(channels{channelToFilter});
block_size=round(1/(2*framerate)*Fs);  %calculate period of 2x signal rate (Nyquist sampling) and how many datapoints that would consist of
thresh=1;  %1 Volt threshold

sig_idx=[];
for i=1:block_size:length(stim_time)-block_size
    [max_val max_location] = max(stim_time(i:i+block_size-1));
    if max_val > thresh
       sig_idx=[sig_idx i+max_location];
    end
end

button = questdlg({['Is ' num2str(numel(sig_idx)) ' the correct number of frames?']},'Frame information','Yes','No','Yes');
if strcmp(button,'No')
    disp('Subtracting last frame...')
    sig_idx1=sig_idx(1:end-1);
end
disp(numel(sig_idx))  %should equal the number of frames
sig_idx1=sig_idx1-1;

figure();
plot(data(1,:));
hold;
plot(sig_idx1,data(1,sig_idx1),'or');

times=convIndex2Time(channels{channelToFilter}, sig_idx1');
% disp(times)
times=ceil(times);
% disp(times)

else 
%     times=[0; channels{channelToFilter}.tim(1:end-1)];  %if CCD exposure signal from PCO pixelfly camera is used with rising edge, frame triggers will be shifted right, towards end of ea exposure. Remove last time and cat 0 to shift frames back left.
%     times = [0; channels{channelToFilter}.tim(1:end-2)];  %if CCD readout signal from PCO pixelfly camera is used with falling edge, there will be +1 frame periods. This will remove the extra readout signal, and shift frame triggers back left by concatenating 0.
    times = channels{channelToFilter}.tim(1:end-1);  %if CCD readout signal from PCO pixelfly camera is used with falling edge, there will be +1 frame periods. This will remove the extra readout signal
    sig_idx1=[];
    disp(numel(times))
end




%{
sig1 = max(stim_time) - (max(stim_time) .* tolerance1); %max signal, with 1% or 5% tolerance
sig_idx1 = find(stim_time>sig1);
sig_pks1 = sig_idx1;

figure();
plot(data(1,:));
hold;
plot(sig_pks1,data(1,sig_pks1),'or');

sig2 = min(stim_time) - (min(stim_time) .* tolerance2); %min signal, with 1% or 5% tolerance
sig_idx2 = find(stim_time<sig2);
sig_pks2 = sig_idx2;

plot(sig_pks2,data(1,sig_pks2),'or');
hold off

times1=sig_pks1;
times2=sig_pks2;
%}


%{


result = zeros(length(sig_idx1),2);
if length(sig_idx1) == length(sig_idx2)
if sig_idx1(1,1) < sig_idx2(1,1)
    
    result(:,1) = sig_idx1';
    result(:,2) = sig_idx2';
else
    result(:,1) = sig_idx2';
    result(:,2) = sig_idx1';
end
else
    error('No. of detected falling and rising edges not equal-- adjust tolerance values!')
end

times=zeros(size(result));
if sum(result,1) ~= 0
times(:,1)=convIndex2Time(channels{channelToFilter}, result(:,1));
times(:,2)=convIndex2Time(channels{channelToFilter}, result(:,2));
disp(times)
tmp=ceil(times);
disp(tmp)
end
%}