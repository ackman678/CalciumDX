function [stn, decpt] = calciumdxdetartifacts(x,nt)
% num=50; %testing
% x=trSign*region.traces(num,:); %testing

% x=mean(ntFilt2,1); %testing

sd=2;  %number of std dev to make threshold for detection
% tr = myfilter(x,2); %perform low pass hann filtering of data with order 2
tr=x;

%----setup defaults for block size and F0 baseline window size-------------
%old def for inmed data (~10Hz in vitro fura2-AM imaging)
% block_size = 3;
% start_baseline = 50;
% end_baseline = 5;


block_size = 3; %delta spacing for spikes and length of window in which to look for local minima
start_baseline = 50;
end_baseline = 5;
%----END setup defaults ---------------------------------------------------

figure;  %testing
ax(1) = subplot(2,1,1); %testing
imagesc(nt) %testing
ax(2) = subplot(2,1,2); %testing
plot(x); %testing
hold on; %testing
linkaxes(ax,'x');

isused = zeros(1,length(tr));
pk = 1;
st = 1;

avgchange = mean(abs(diff(x)));
stdchange = std(abs(diff(x)));
baseavgs=[];

%--try to get initial st (spike time) and pk (peak time) estimates-----------------------
for c = 3:block_size:length(tr)-block_size  %loop through block sizes
    [dummy min_location] = min(tr(c:c+block_size-1)); %fetch local minima (our possible signals)
    min_location = min_location+c-1;
    
    %setup baseblock and blocksize
    if min_location - pk(end) < 2*block_size  %if next min_location is inside twice the block size from last pk
        baseblock = tr(pk(end):min_location); %baseblock will be from last peak to current min_location
        baseavg = max(baseblock);
    else %if next min_location is far enough from last pk
        baseblock = tr(max([1 min_location-start_baseline]):max([1 min_location-end_baseline]));
        %         baseavg = mean(baseblock(baseblock>=median(baseblock))); %problem is here in baseavgg. Weighted mean for base avg. Take the mean of those values above the median. Baseavg get shifted if there are any large positive peak in the baseblock.
        baseavg = median(baseblock); %much better. Now not susceptible to positive dF influences. This should be okay, especially for already weighted/smoothed data (like detection of zartifact transients on averaged tr data)
    end
    
    %     avgchange = mean(abs(x(2:end)-x(1:end-1)));  %same as mean(abs(diff(x)))  Doesn't need to be inside this for loop moved above
    %     stdchange = std(abs(x(2:end)-x(1:end-1))); %same as std(abs(diff(x)))  Doesn't need to be inside this for loop moved above
    
    %no. of sd for detection set here
    if tr(min_location)-baseavg < -(avgchange+sd*stdchange) %problem is in baseavg
        ps = find(tr>baseavg); %look here for adjustment %problem is in baseavg
        ps = ps(ps < min_location); %look here for adjustment
        if ~isempty(ps)
            ps = ps(end); %single point set as spike onset first happens here
            if ps < pk(end)  %this is used to make an adjustment if the current spike onset was found to overlap with the signal for the last spike. Pretty necessary for this algorithm.
                [dummy psn] = max(tr(pk(end)+1:min_location)); %set psn as a new possible spike onset location.
                psn = psn+pk(end); %set the new possible spike onset to the appropriate index number.
                if tr(min_location) < tr(psn)
                    ps = psn;  %keep this spike onset
                else
                    ps = psn;
                    isused(ps) = 1;
                end
            end
            if isused(ps) == 0  %if this onset hasn't been used yet, add it the spike time and correponding peak times to the st and pk vectors.
                %                 plot([min_location ps],[tr(min_location) tr(ps)],'-r'); %testing
                pk = [pk min_location]; %local minima accepted as signal peak at this baseblock location
                st = [st ps]; %corresponding spike time from ps calculation set here.
                baseavgs=[baseavgs baseavg];
            end
        end
    end
end
%}

% if pk is more than one data point keep pk and st, otherwise set to empty because first pk point was set to '1' as a placeholder.
if length(pk>1)
    pk = pk(2:end);
    st = st(2:end);
else
    pk = [];
    st = [];
end


%setup stn(c) and pkn(c) based on st(c) and pk(c) if st(c)-pkn(end) is less than block_size threshold. Or if pk is empty, set pkn and stn to empty.
if~isempty(pk)
    pkn = pk(1);
    stn = st(1);
    baseavgs2=baseavgs(1);
    for c = 2:length(pk)
        if st(c)-pkn(end) <= block_size %if time between a peak and subsequent spike is too little (less than block size) then merge them.
            pkn(end) = pk(c);
            %             baseavgs=[baseavgs(1:c-1) baseavgs(c+1:end)];
        else
            stn(end+1) = st(c);
            pkn(end+1) = pk(c);
            baseavgs2=[baseavgs2 baseavgs(c)];
        end
    end
else
    pkn = [];
    stn = [];
end


%Old stn(c) refinement-----------------------------------------------------
% for c = 1:length(pkn)
%     th = tr(stn(c))-(tr(stn(c))-tr(pkn(c)))/50;
%     f = find(tr>th);
%     f = f(find(f<pkn(c)));
%     if ~isempty(f)
%         stn(c) = f(end);
%     end
% end



%New stn(c) refinement------------------------------------------------------
tr = x; %use the original non-filtered data
for c = 1:length(pkn)
    startpt = max([1 stn(c)-end_baseline]);  %just in case the current spike time is near beginning of recording. Otherwise end_baseline variable is the lookback index location
    endpt = find(tr>tr(stn(c))-(tr(stn(c))-tr(pkn(c)))/5);  %find all trace data greater than a 20% fractional change from the overall signal dF for the current minima peak.
    endpt = endpt(endpt<pkn(c)); %set endpts indices to those that occur before the current peak index
    endpt = endpt(end)+1; %set endpt to this last endpt
    dfblock = tr(startpt:endpt); %this is now a block of time surrounding the current stn(c) spike onset
    %the following lines iteratively refine down to where the spike onset should be located.
    f = find(dfblock(2:end-1)>dfblock(3:end) & dfblock(2:end-1)>dfblock(1:end-2))+startpt;
%     disp(num2str(c)) %for troubleshooting
%     disp(num2str(f)) %for troubleshooting
    if ~isempty(f)
        if length(tr)-f(end)-3>0  % paolo added
            f = f(end);
            for fx = 1:2
                if tr(f+fx+1)-tr(f+fx) < 4*(tr(f+fx)-tr(f))
                    f = f+fx;
                end
            end
            stn(c) = f;
        end % paolo added
        stn(c) = f(end);
%         disp(num2str(stn(c)))  %for troubleshooting
    end
    %addded by JBA 2010-10-23
    if tr(stn(c)) > baseavgs2(c)+3*stdchange
        indRange=stn(c):pkn(c);
        ind=find(tr(indRange)>baseavgs2(c));
        %         ind=find(tr(indRange)>median(baseavgs));
        %         ind=find(indRange<mean(baseavgs));
        ind=ind+(stn(c)-1);
        if ind<1
            ind=1;
        end
        stn(c)=ind(end);
%         disp(num2str(stn(c))) %for troubleshooting
    end
end

%}


%refine which spike times (stn) we keep--
%keeping only those peaks that are at least 20% in amplitude of the largest ampl transient in the recording. For deleting small ampl transients (like from a second A.P. or SPA cell fluctuations) that are riding on top of another, much larger transient
%{
f = find(tr(stn)-tr(pkn) > 0.2*max(tr(stn)-tr(pkn)));
stn = stn(f);
pkn = pkn(f);
%}

%fetch offset points based on where stn(c) is located-- what is the '15' constant? Spike onset + 15frames? Yes. Or 1.5 sec at 10fps, the framerate of hardware originally used with this code.
decpt = [];
for c = 1:length(stn)
    [mn i] = min(tr(stn(c):min([stn(c)+15 size(x,2)])));  %find the minimum point index between stn(c) and 15 frames later (or end of recording). Should be signal peak.
    i = i+stn(c)-1; %set signal peak index i to the appropriate frame number.
    %     th = (tr(stn(c))+tr(i))/2; %offset will be half of dF. Could set this threshold value to something different, like getting entire transient period (set within 1SD of the total baseline signal)
    th = (tr(stn(c))+tr(i))/4; %set offset to within one-quarter of dF.
    f = find(tr>th); %fetch all trace indices greater than the threshold
    f = f(f>i); %limit it to those indices occuring after the signal peak index i.
    if isempty(f)
        decpt(c) = size(x,2);
    else
        decpt(c) = f(1); %set the offset to the first point occuring at half of peak value after the signal peak has occured
    end
    if c < length(stn) && decpt(c) >= stn(c+1) %if c is not near end of recording and the current index is overlapping with the next spike onset index then make the current offset set to one frame before the next onset.
        decpt(c) = stn(c+1) - 1;
    end
    %Could add a loop here that merges transients if the time between offsets and subsequent onsets is less than blocksize
end

%Testing----------------------------------------------
for c = 1:length(pkn)
    plot([stn(c) pkn(c)],[tr(stn(c)) tr(pkn(c))],'-k');
    plot(stn(c),tr(stn(c)),'+k');
    plot(decpt(c),tr(decpt(c)),'*k');
end

% for c = 1:length(pk)
%     plot([st(c) pk(c)],[tr(st(c)) tr(pk(c))],'-r');
% end
%End Testing------------------------------------------
%}