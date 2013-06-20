%corrWavesBilateral.m
%script to measure retinal wave correlation between regions (hemispheres), based on temporal-spatial-direction properties of the observed wave distribution. Corr measusred using Baysian inference/Monte Carlo simulations
%input= filelist
%James B. Ackman 2011-12-22

function results = corrWavesBilateralSimTcorr(filelist,region,numres)
%results=corrWavesBilateralSimTcorr({filename},region,1000);  %for a single file
%results=corrWavesBilateralSimTcorr(filelist,region,1000);  %for a list of files
%==(3)==  Randomize wave interval spike times=============================================

makefigs = 1;

%LOAD FILE
if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end
%results={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency'};
fnms=filelist(:,1);
for fIdx=1:numel(fnms)
    
    if loadfile > 0
        myMATfile=load(fnms{fIdx});
        region=myMATfile.region;
    end
    [pathstr, name, ext] = fileparts(fnms{fIdx});
    rowinfo = {[name ext]};  %2011-07-11 jba
    %     rowinfo = filelist(j,:);
    sprintf(fnms{fIdx})
    h = waitbar(fIdx/numel(fnms));
    
    %----SETUP PARAMETERS----
    %myParams.waveinterval_mu = 65.15658; %calculated from the observed waveisiOns dataset in R
    %myParams.waveinterval_sigma = 46.3379; %calculated from the observed waveisiOns dataset in R
    param.numres = numres;
    param.deltaspacing=10; %10secs in between waves was used as a refractory time in detecting waves in calciumdxDetectWaves.m
    param.maxt=size(region.traces,2)*region.timeres; %600s in between waves was used as a refractory time in detecting waves in calciumdxDetectWaves.m
    param.locationIndices = [2 3];  %the length should be the same as the value of nRegions if both hemisphere spike times want to be reshuffled
    param.nRegions = [1];  %should match the no. of unique locationIndices if both hemisphere spike times want to be reshuffled, otherwise must be less than the length of param.locationIndices, like a value of 1 to reshuffle the second hemipshere.
    %	param.regionnames(1) = 'SC.R';
    %	param.regionnames(2) = 'SC.L';
    param.movieLength = size(region.traces,2); %should be set by the frame length of movie currently being analysed
    param.nwaves = [length(region.wavedata{param.locationIndices(1)}.waveonsets) length(region.wavedata{param.locationIndices(2)}.waveonsets)];  %should be set by the observed no. of waves in the movie currently being analysed
    param.framePeriod = region.timeres;
    param.sig = 15;  %sigma for gaussian smoothing.  No. of +/- frames for the 37% mark in the PDF for the gaussian function.  i.e. for sig = 1, resulting spike matrix will be [0.0001 0.0183 0.3679 1.0000 0.3679 0.0183 0.0001] surrounding the spike onset.
    param.p_val = 0.05;
    param.constrainNwaves = 0;
    param.fnmbase = [name '_xcorrn_' datestr(now,'yyyymmdd-HHMMSS')];  %base filename for saving graphics output
    %GET OBSERVED correlation matrix-----------------------------
    s = setupWaveTimeMatrix(region,param);
    s_thick = gauss_events(s,param.sig);
    %figure; imagesc(s_thick); title('obs')
    obs.corr = s_thick*s_thick';
    obs.corr_cont = reshape(obs.corr,1,numel(obs.corr)); %reshape adjacency matrix into a vector
    
    %DO SIMULATIONS on randomized data====================================================
   
    npairs = zeros(2,param.numres);
    countn = zeros(2,param.numres);  %setup the results histogram array
    ncorrcoef = zeros(2,param.numres);
    
    
    hbar = waitbar(0,'Please wait...');
    if makefigs > 0
    scrsz = get(0,'ScreenSize');
    fig1=figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(3)/2]);
    end
    for t = 1:param.numres
        if rem(t,10) == 0
            waitbar(t/numres,hbar);
        end
        rand_s = setupRandomWaveTimeMatrix(param,s);
        rand_sthick = gauss_events(rand_s,param.sig);
        %figure; imagesc(rand_sthick); title('res')
        res.corr = rand_sthick*rand_sthick';
        
    	[r,p] = corrcoef(rand_sthick');
        ncorrcoef(1,t) = p(2,1);
        ncorrcoef(2,t) = p(2,1);
        
        tmpIND = find(rand_sthick<0.01);
        rand_sthick(tmpIND) = 0;
        
        if makefigs > 0
            set(0,'CurrentFigure',fig1)
            imagesc(rand_sthick); colormap(jet); %colorbar
            M1(t) = getframe;
        end
        
        
        
        %--PUT HERE-- get obs.corr spike indices-------------------------------
        i=2;
        j=1;
        %find the correlated spike indices for cell i, cell j
        cellIND = rand_sthick(i,:) .* rand_sthick(j,:);  %element wise multiplication to intersect the spike times for the two corr cells
        [muhat, sigmahat] = normfit(cellIND);
        level = muhat+(2*sigmahat);
        spkIND = find(cellIND>level);
        if ~isempty(spkIND)
            dfSpk = find(diff(spkIND) > 1);
            if ~isempty(dfSpk)
                spkIND = setdiff(spkIND,spkIND+1)   %will output the spk indices that are separated (the no. of spks)
            else
                spkIND = round(mean(spkIND))   %if there is just one gauss spk found, then make equal to middle of spk
            end
            
            spkINDlength = length(spkIND);
            npairs(1,t) = spkINDlength;
            npairs(2,t) = spkINDlength;
            countn(1,t) = length(find(rand_s(1,:)));
            countn(2,t) = length(find(rand_s(2,:)));
            
        else
            npairs(1,t) = 0;
            npairs(2,t) = 0;
            countn(1,t) = length(find(rand_s(1,:)));
            countn(2,t) = length(find(rand_s(2,:)));
        end
        
        %=================================================================================
    end  %end of for t loop
    close(hbar)

    %-------------------------------------------
    %data=results;
    results(fIdx).filename = fnms{fIdx};
    results(fIdx).countn = countn;
    results(fIdx).date=datestr(now,'yyyymmdd-HHMMSS');
    results(fIdx).npairs = npairs;
    results(fIdx).ncorrcoef= ncorrcoef;
    results(fIdx).obs = obs;
    results(fIdx).obs.s = s;
    results(fIdx).obs.s_thick = s_thick;
    results(fIdx).res = res;
    results(fIdx).res.rand_sthick = rand_sthick;
    results(fIdx).param = param;
%     figure; imagesc(results(1).obs.s_thick); title(['obs' ', GaussSig=' num2str(param.sig)]); zoom xon  %TESTING
    figure; imagesc(results(1).res.rand_sthick); title(['res at t=1000' ', GaussSig=' num2str(param.sig)]); zoom xon  %TESTING
%     myCDF(npairs,'percent pairs'); title(['npairs over time' ', GaussSig=' num2str(param.sig)])  %TESTING
    save([param.fnmbase '.mat'],'results')    
    
    %---Printout information---
    disp('No. of waves overlapping')
    sum(results.npairs,2)
    
    disp('No. of random waves')
    sum(results.countn,2)
    
    disp('fraction of waves overlapping')
    sum(sum(results.npairs,2))/sum(sum(results.countn,2))

    [m,n]=find(results.ncorrcoef < 0.05);
    disp(numel(m))
    
        if makefigs > 0
            vidObj = VideoWriter([param.fnmbase '-rand_sthick' '.avi'])
            open(vidObj)
            for i =1:numel(M1)
            writeVideo(vidObj,M1(i))
            end
            close(vidObj)
        end
    
end
close(h)
end

%=========================================================================================
function s = setupWaveTimeMatrix(region,param)
s = zeros(length(param.locationIndices),param.movieLength);
for i = 1:length(param.locationIndices)
    waveonsets = region.wavedata{param.locationIndices(i)}.waveonsets;
    s(i,waveonsets) = 1;
end
%s_thick = gauss_events(s,param.sig);
%figure; imagesc(s_thick); title('obs')
end

%=========================================================================================
function s = setupRandomWaveTimeMatrix(param,s_obs)
%--set up random spike matrix dataset-------
s = zeros(size(s_obs));
for i = param.nRegions
    if length(param.nRegions) < length(param.locationIndices)
        s(i,:) = s_obs(i,:);  %set the first hemisphere to the obs hemsiphere wave times
        i = i+1;   %if we are only reshuffling one hemisphere, then it will be the second ROI.
    end
    randWaveIntervals = getRandomWaveIntervals(param.nwaves(i), [], [], param.maxt,[],param.constrainNwaves);
    randWaveTimes=cumsum(randWaveIntervals);
    %disp(randWaveTimes)
    %disp(sum(randWaveIntervals))
    if ceil(randWaveTimes(end)/param.framePeriod) > param.movieLength && param.constrainNwaves > 0
        randWaveTimes = randWaveTimes-param.deltaspacing;
        s(i,ceil(randWaveTimes / param.framePeriod)) = 1;
    elseif param.constrainNwaves > 0
        s(i,ceil(randWaveTimes / param.framePeriod)) = 1;
    else
        wvIdx = find(randWaveTimes < param.maxt);
        s(i,ceil(randWaveTimes(wvIdx) / param.framePeriod)) = 1;
    end
    %error(num2str(randWaveTimes(end)))
end
%rand_sthick = gauss_events(s,param.sig);
%figure; imagesc(rand_sthick); title('res')
end


