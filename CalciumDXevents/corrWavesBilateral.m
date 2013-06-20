%corrWavesBilateral.m
%script to measure retinal wave correlation between regions (hemispheres), based on temporal-spatial-direction properties of the observed wave distribution. Corr measusred using Baysian inference/Monte Carlo simulations
%input= filelist
%James B. Ackman 2011-11-16
%==(1)==  Get distribution of interwave intervals (onset-onset) waveonset.fr in database==
%in R  http://www.r-project.org
%#=====2011-11-16---New interwave intervals based on onset-onset distances-------
%quartz(); qplot(waveisi.s,data=dwaves2,geom="histogram") + xlab("Waveisi.s, off-on")  #existing waveisi.s distribution
%regionsToMatch <- as.character(with(dwaves2,filename:region.name))  #filename by hemisphere
%filelistMatch <- with(dwaves2,unique(regionsToMatch))
%dwaves2$waveisiOns.s <-  dwaves2$waveisi.s
%for (i in 1:length(filelistMatch)) {
%idx <- with(dwaves2,which(regionsToMatch == filelistMatch[i]))
%frDiff <- with(dwaves2,diff(waveonset.fr[idx]))
%timeDiff <- with(dwaves2,frDiff*framePeriod[idx[1]])
%dwaves2$waveisiOns.s[idx[-length(idx)]] <- timeDiff
%}
%m <- dwaves2; quartz(); qplot(waveisiOns.s,data=m,geom="histogram") + xlab("Waveisi.s, on-on")  #new. 
%quartz(); qplot(log(waveisiOns.s),data=m,geom="histogram") + xlab("Waveisi.s, on-on")  #looks pretty normal in the log transformation
%
%==(2)==  Fit waveisiOns.s to normal distribution, obtain mu and sigma====================
%the log transformed waveisiOns data set is gaussian distributed
%shapiro.test(log(m$waveisiOns.s)) #data is normally distributed
%qqnorm(log(m$waveisiOns.s)); qqline(log(m$waveisiOns.s),col=2);  #data is normally distributed
%mu <- mean(m$waveisiOns.s,na.rm=TRUE)  #gives mean of 70.1123 s
%sigma <- sd(m$waveisiOns.s,na.rm=TRUE) #gives sd of 52.56837 s.
%deltaspacing=10; #10secs in between waves
%randWaveIntervals <- rnorm(1000,mean=mu,sd=sigma)
%answer <- which(randWaveIntervals < 10)
%while (length(answer) > 0) {
%	newInt <- abs(rnorm(length(answer),mean=mu,sd=sigma))
%	randWaveIntervals[answer] <- newInt
%	answer <- which(randWaveIntervals < 10)
%	}
%quartz(); qplot(randWaveIntervals,data=m,geom="histogram") + xlab("RandN wave intervals"); ggsave(file="hist_randNwaveInts.pdf")

%Fit waveisisOns.s to gamma distribution
%x <- na.omit(dwaves2$waveisiOns.s)
%fitdistr(x,"normal")
%     # mean         sd    
%  # 70.112296   52.539219 
% # ( 1.749365) ( 1.236988)
%fitdistr(x,"gamma")
%      # shape         rate    
%  # 2.355830701   0.033592751 <-----------------scale = 1/rate, matlab requires scale, B
% # (0.103783129) (0.001647333)
%       # shape         rate    #for just CCD.img waves
%  # 2.517740215   0.038645565 
% # (0.117240497) (0.001989708)

%Ashape = 2.517740215;  %smaller shape values shift distribution left    %TESTING
%Bscale = 1/0.038645565;   %TESTING
%R = 10 + gamrnd(Ashape,Bscale,1000, 1) % A is shape parameter, B is scale parameter   %TESTING
%figure; hist(R,20)  %TESTING
%randWaveIntervals = getRandomWaveIntervals(1000, myParams.waveinterval_mu, myParams.waveinterval_sigma, myParams.deltaspacing); %Testing Distribution
%hist(randWaveIntervals,20)  %TESTING

function results = corrWavesBilateral(filelist,region,numres,wavesDBfnm)
% filelist = readtext('testfiles_20120103.txt',' ');
%results=corrWavesBilateral({filename},region,1000,'wavesDB.mat');  %for a single file
%results=corrWavesBilateral(filelist,[],1000,'wavesDB.mat');  %for a list of files
%==(3)==  Randomize wave interval spike times=============================================

%LOAD FILE
if nargin< 2 || isempty(region); region = []; end
if isempty(region)
    loadfile = 1;
else
    loadfile = 0;
end
%results={'matlab.filename' 'region.name' 'roi.number' 'nrois' 'roi.height.px' 'roi.width.px' 'xloca.px' 'yloca.px' 'xloca.norm' 'yloca.norm'  'stimulus.desc' 'nstimuli' 'activeResponseFraction.region' 'responseFreq.norm' 'absoluteFiringFreq.Hz' 'meanLatency'};
fnms=filelist(:,1);
fnms2=filelist(:,2);
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
	param.nRegions = 1;  %should match the no. of unique locationIndices if both hemisphere spike times want to be reshuffled, otherwise must be less than the length of param.locationIndices, like a value of 1 to reshuffle the second hemipshere.
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
	
	%--PUT HERE-- get obs.corr spike indices-------------------------------
	i=2;
	j=1;
	%find the correlated spike indices for cell i, cell j
	cellIND = s_thick(i,:) .* s_thick(j,:);  %element wise multiplication to intersect the spike times for the two corr cells
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
	
		%find the obs overlapping wave indices-----------------------
		for locationIndex = 1:length(param.locationIndices)
			tmp = [];
			for spkINDlength = 1:length(spkIND)
				[dummy, wvIND] = min(abs(region.wavedata{param.locationIndices(locationIndex)}.waveonsets - repmat(spkIND(spkINDlength),1,length(region.wavedata{param.locationIndices(locationIndex)}.waveonsets))));
				tmp = [tmp wvIND];
			end
			obs.waves(locationIndex).waveIND = tmp;
			obs.waves(locationIndex).regionname = region.wavedata{param.locationIndices(locationIndex)}.name;
		end
	else
		obs.waves(locationIndex).waveIND = [];
		obs.waves(locationIndex).regionname = region.wavedata{param.locationIndices(locationIndex)}.name;
	end

	%get wavesDB lookup table-----------------------------------			
%	wavesDBfnm = 'wavesDB.mat';
	matObj = matfile(wavesDBfnm);
	szDB = size(matObj,'filename',2);
	waveDBidx = [1:szDB]';
	tmp = matObj.filename;
	filename = {tmp.filename}';
	tmp = matObj.regionname;
	regionname = {tmp.regionname}';
	tmp = matObj.wavenumber;
	wavenumber = [tmp.wavenumber]';
%	for i = 1:szDB
%		disp([num2str(waveDBidx(i)) ' ' filename{i} ' ' num2str(regionname{i}) ' ' num2str(wavenumber{i})])
%	end

	%--PUT HERE-- get obs.spatialcorr metrics---------------------------------------------	
%	testFilename = '/Volumes/FIRMTEK/110323i/110323_08.tif';  %for TESTING on 2011-12-12
	testFilename = fnms2{fIdx};  %use this
	disp(testFilename)
	%temporary blacklist for waves in wavesDB not to include in the random reshuffling script below--
	blacklist = {'/Volumes/FIRMTEK/110323i/110323_01.tif' '/Volumes/FIRMTEK/110323i/110323_03.tif' '/Volumes/FIRMTEK/110323i/110323_04.tif'}
	blacklistIdx = [];
	for blklistCounter = 1:numel(blacklist)
		blacklistIdx = [blacklistIdx; find(strcmp(filename,blacklist{blklistCounter}))];
	end
	disp(obs.waves(1).waveIND)
	obs.Mahalmetrics =[];
	obs.Euclmetrics = [];
	for wvIND = 1:length(obs.waves(1).waveIND)
		disp(wvIND)		
		%use wavesDB lookup table to get waveDB indices 
		wave1number = find(strcmp(filename,testFilename) & strcmp(regionname,obs.waves(1).regionname) & wavenumber == obs.waves(1).waveIND(wvIND));
		wave2number = find(strcmp(filename,testFilename) & strcmp(regionname,obs.waves(2).regionname) & wavenumber == obs.waves(2).waveIND(wvIND));
		waveCorrRefno = wave2number;	%use wave in hemi 2 as reference
		wave1number
		wave2number

		
		xcorrResults = testXcorrDistMetric(wavesDBfnm,wave1number,wave1number);   %autocorr reference
		obs.data(wave1number).xcorrResults = xcorrResults;
		xcorrResults = testXcorrDistMetric(wavesDBfnm,wave1number,wave2number);   %obs corr 
		obs.data(wave2number).xcorrResults = xcorrResults;
		
		obs.wavesDBind(wvIND).pair = [wave1number wave2number];
		obs.wavesIND(wvIND).pair = [obs.waves(1).waveIND(wvIND) obs.waves(2).waveIND(wvIND)];
		obs.filename(wvIND).filename = testFilename;

		[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(obs.data,wave1number,wave2number,waveCorrRefno);				
		obs.data(waveCorrRefno).DistMahal = DistMahal;
		obs.data(waveCorrRefno).meanDistMahal = meanDistMahal;
		obs.data(waveCorrRefno).AllSqDistFromCenter = AllSqDistFromCenter;
		obs.data(waveCorrRefno).meanEuclDist = meanEuclDist;
		
		[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(obs.data,wave1number,wave1number,waveCorrRefno);   %get the orig ref for comparison.				
		obs.data(wave1number).DistMahal = DistMahal;
		obs.data(wave1number).meanDistMahal = meanDistMahal;
		obs.data(wave1number).AllSqDistFromCenter = AllSqDistFromCenter;
		obs.data(wave1number).meanEuclDist = meanEuclDist;
		
		obs.Mahalmetrics =[obs.Mahalmetrics abs(obs.data(waveCorrRefno).meanDistMahal - obs.data(wave1number).meanDistMahal)];
		obs.Euclmetrics = [obs.Euclmetrics obs.data(waveCorrRefno).meanEuclDist];
	end
	%--------------END get obs.corr spks--------------------------------------------------






	
	%DO SIMULATIONS on randomized data====================================================
	corr_res = zeros(fix(param.p_val*param.numres)+1,numel(obs.corr));  %this line uses the pval. Less rows from smaller p_vals (probability of event in no. of simulations) or less num_reshuffles.
	npairs = zeros(1,param.numres);
	countn = zeros(size(obs.corr));  %setup the results histogram array
	hbar = waitbar(0,'Please wait...');
	
	res.DistMetricsCountn = zeros(1,length(obs.waves(1).waveIND));
	res.MahalmetricsCountn = zeros(1,length(obs.waves(1).waveIND)); 
	res.EuclmetricsCountn = zeros(1,length(obs.waves(1).waveIND));
	
	for t = 1:param.numres
		if rem(t,10) == 0
			waitbar(t/numres,hbar); 
		end 
		rand_s = setupRandomWaveTimeMatrix(param,s);
		rand_sthick = gauss_events(rand_s,param.sig);
		%figure; imagesc(rand_sthick); title('res')
		res.corr = rand_sthick*rand_sthick';
	%	[r,p] = corrcoef(rand_sthick');
		%find temporally correlated cell pairs---------------------------------
		countn = countn + (res.corr>=obs.corr);  %from corr_pairs_spatialdir.m   %add to counter...
		
		corr_res(1,:) = reshape(res.corr,1,numel(res.corr)); %add the reshuffled adjacency matrix corr values to the first row of probability array.
	%         pairs = (corrs_cont>corrs_res(1,:));  %jba
		corr_res = sort(corr_res); %sort columns of probablity matrix in order of ascending corr value
		pairs = (obs.corr_cont>corr_res(2,:));  %jba
		pairs = reshape(pairs,size(rand_s,1),size(rand_s,1));
		disp(pairs)
		[i j] = find(pairs==1);
		f = find(j>i);
		pairs = [i(f) j(f)];
		npairs(t) = 100*size(pairs,1)/(size(s,1)*((size(s,1)-1))/2);
				
		
		%==========BEGIN spatial corr loop================================================
		res.pairs = res.corr>=obs.corr;  %from corr_pairs_spatialdir.m
		%If temporally corr cell pairs are present then find spatial and direction 
		%correlated spike times and set up results arrays----------------------
		[iPairs jPairs] = find(tril(res.pairs,-1));  %find non zero counts of pairs in the newly acquired adjacency matrix from the reshuffled spike trains
		if ~isempty(iPairs)
			for c=1:size(iPairs,1)
				%get the corr cell indices i,j
				%         size([i j],1);
				i = iPairs(c);
				j = jPairs(c);
				%find the correlated spike indices for cell i, cell j
				cellIND = rand_sthick(i,:) .* rand_sthick(j,:);  %element wise multiplication to intersect the spike times for the two corr cells
				%             level = graythresh(cellIND);  %separate the histrogram of corr values into two populations, especially useful if gauss events was used for spike smoothing and there is loads of near-zero values.
				[muhat, sigmahat] = normfit(cellIND);
				level = muhat+(2*sigmahat);
				spkIND = find(cellIND>level);
				
				dfSpk = find(diff(spkIND) > 1);
				if ~isempty(dfSpk)
					spkIND = setdiff(spkIND,spkIND+1);
				end
				
				%--PUT HERE-- get res.corr spike indices-----------
				for locationIndex = 1:length(param.locationIndices)
					tmp = [];
					for spkINDlength = 1:length(spkIND)
						[dummy, wvIND] = min(abs(region.wavedata{param.locationIndices(locationIndex)}.waveonsets - repmat(spkIND(spkINDlength),1,length(region.wavedata{param.locationIndices(locationIndex)}.waveonsets))));
						tmp = [tmp wvIND];
					end
					res.waves(locationIndex).waveIND = tmp;
					res.waves(locationIndex).regionname = region.wavedata{param.locationIndices(locationIndex)}.name;
				end	
				%---END get res.corr spike indices-----------------
				
				%DO Random 2D xcorr/distance metric for spatial wave pattern similarity from wavesDB.   %<-------------------------------**
				%BASE the NO. of waves to try corr metric on no. of spkIND above   %<-------------------------------**
			%
				%--PUT HERE-- get res.spatialcorr metrics---------------------------------------------	
				Mahalmetrics =[];
				Euclmetrics = [];
				for wvIND = 1:length(obs.waves(1).waveIND)  %<--** %first test repeat for each obs overlapping wave to add to counter... Alternatively, change to no. of reshuffled overlapping waves, use randn to generate which of the obs. overlapping wave indices to make the xcoor comparison below for. Make counters for total no. of sign corr 
					%use wavesDB lookup table to get waveDB indices 
					wave1number = find(strcmp(filename,testFilename) & strcmp(regionname,obs.waves(1).regionname) & wavenumber == obs.waves(1).waveIND(wvIND));
					% wave2number = find(strcmp(filename,testFilename) & strcmp(regionname,obs.waves(2).regionname) & wavenumber == obs.waves(2).waveIND(wvIND));
					% waveCorrRefno = obs.waves(2).waveIND(wvIND);	%use wave in hemi 2 as reference
					
					randwave2number = randi(szDB,1,1);	
					while (randwave2number == wave2number) || (ismember(randwave2number,blacklistIdx))
						randwave2number = randi(szDB,1,1);   %obtain a rand wave no. from the wavesDB
					end
					waveCorrRefno = randwave2number;
					disp(wvIND)
					disp(['random wave2number: ' num2str(randwave2number)])
					
					% xcorrResults = testXcorrDistMetric('wavesDB.mat',wave1number,wave1number);   %autocorr reference, don't need to run this a second time.
					res.data(wave1number).xcorrResults = obs.data(wave1number).xcorrResults;   %since we already have the autocorr baseline stored here...
					xcorrResults = testXcorrDistMetric(wavesDBfnm,wave1number,randwave2number);   %obs corr 
					res.data(randwave2number).xcorrResults = xcorrResults;
					
					res.wavesDBind(wvIND).pair = [wave1number randwave2number];
					%res.wavesIND(wvIND).pair = [obs.waves(1).waveIND(wvIND) obs.waves(2).waveIND(wvIND)];   %only makes sense if we change to the Alternative wvIND fetching above...
					res.filename(wvIND).filename = testFilename;
			
					[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(res.data,wave1number,randwave2number,waveCorrRefno);				
					res.data(waveCorrRefno).DistMahal = DistMahal;
					res.data(waveCorrRefno).meanDistMahal = meanDistMahal;
					res.data(waveCorrRefno).AllSqDistFromCenter = AllSqDistFromCenter;
					res.data(waveCorrRefno).meanEuclDist = meanEuclDist;
					
					[DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist] = getSpatialCorrDistanceMetrics(res.data,wave1number,wave1number,waveCorrRefno);   %get the orig ref for comparison.				
					res.data(wave1number).DistMahal = DistMahal;
					res.data(wave1number).meanDistMahal = meanDistMahal;
					res.data(wave1number).AllSqDistFromCenter = AllSqDistFromCenter;
					res.data(wave1number).meanEuclDist = meanEuclDist;
					
					mDist = abs(res.data(waveCorrRefno).meanDistMahal - res.data(wave1number).meanDistMahal);
					Mahalmetrics =[Mahalmetrics mDist];
					Euclmetrics = [Euclmetrics res.data(waveCorrRefno).meanEuclDist];
					
					if (mDist <= obs.Mahalmetrics(wvIND)) || (res.data(waveCorrRefno).meanEuclDist <= obs.Euclmetrics(wvIND))   %here is the test statistic. Add to counter and calculate pval based on param.numres below.
						res.DistMetricsCountn(wvIND) = res.DistMetricsCountn(wvIND) + 1;					
						if (mDist <= obs.Mahalmetrics(wvIND))
						res.MahalmetricsCountn(wvIND) = res.MahalmetricsCountn(wvIND) + 1;
						elseif (res.data(waveCorrRefno).meanEuclDist <= obs.Euclmetrics(wvIND))
						res.EuclmetricsCountn(wvIND) = res.EuclmetricsCountn(wvIND) + 1;
						end
					end
					
				end   %end of wvIND for loop
				%}
			end   %end of size iPairs loop
		res.Mahalmetrics(t).Mahalmetrics = Mahalmetrics;  %save the measured value for each simulation loop
		res.Euclmetrics(t).Euclmetrics = Euclmetrics;  %save the measured value for each simulation loop
		end  %end of if spatial loop		
		%=================================================================================
			
    end  %end of for t loop
    close(hbar)
    fprintf('\n');
    if ~isempty(pairs)
    disp('significant temporal pairs at locations: ')   %Are the hemispheres sign t corr?
    disp(pairs)
    else
    disp('no significant temporal pairs')
    end
    p_val = countn./param.numres;
    tmpIdx = find(tril(ones(size(p_val)),-1));
    disp('temporal p = ')
    disp(num2str(p_val(tmpIdx)'))
    
    disp('Combined Spatial corr counts:')
    disp(num2str(res.DistMetricsCountn))
    disp('spatial pvals:')
    DistMetricsPvals = res.DistMetricsCountn./param.numres; 
    disp(num2str(DistMetricsPvals))
    disp('Mahal Spatial corr counts:')
    disp(num2str(res.MahalmetricsCountn))
    disp('Eucl Spatial corr counts:')
    disp(num2str(res.EuclmetricsCountn))
    
    disp('Corresponding wavesDB indices:')
    disp(vertcat(obs.wavesDBind(1:length(obs.waves(1).waveIND)).pair))
    disp('Corresponding file wave numbers:')    
    disp(vertcat(obs.wavesIND(:).pair))
	%Tally up and display counter information and pvals for:
	%How many waves in movie have t corr?
	%Do the hemispheres have sign TDL corr?
	%Which/how many waves in movie have sign TDL corr?
 
    %-------------------------------------------
	%data=results;
	results(1).filename = fnms{fIdx};
	results(1).countn = countn;
	results(1).temporal_p_val = p_val(tmpIdx)';
	results(1).date=datestr(now,'yyyymmdd-HHMMSS');
	results(1).npairs = npairs;
	results(1).pairs = pairs;
	results(1).obs = obs;
	results(1).obs.s = s;
	results(1).obs.s_thick = s_thick;
	results(1).res = res;
	results(1).res.rand_sthick = rand_sthick;
	results(1).param = param;
%	figure; imagesc(results(1).obs.s_thick); title(['obs' ', GaussSig=' num2str(param.sig)]); zoom xon  %TESTING
%	figure; imagesc(results(1).res.rand_sthick); title(['res at t=1000' ', GaussSig=' num2str(param.sig)]); zoom xon  %TESTING
%	myCDF(npairs,'percent pairs'); title(['npairs over time' ', GaussSig=' num2str(param.sig)])  %TESTING
	save([param.fnmbase '.mat'],'results')
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
for i = 1:param.nRegions
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

%=========================================================================================
function [DistMahal, meanDistMahal, AllSqDistFromCenter, meanEuclDist]  = getSpatialCorrDistanceMetrics(data,wave1number,wave2number,waveCorrRefno)
	%--PUT HERE-- test the spatial corr metric--
%	allmeanDistMahal  = [];
%	allmeanEuclDist = [];
	if ~isempty(data(waveCorrRefno).xcorrResults.dXY)
	nRows = size(data(waveCorrRefno).xcorrResults.dXY,1);  %for using observed reference wave as reference
	randJitter = rand(nRows,2);  %for using observed reference wave as reference
	dXY1 = data(waveCorrRefno).xcorrResults.dXY + randJitter;  %for using observed reference wave as reference
	end
%	for i = wave2numbers
%		disp(num2str(i))
	%     i=27
		i = wave2number;
		regpropInd = [1:10];  %in case you want to limit the no. of properties input to the calculation from regionprops
		textureInd = [1:2 4:6];
	%     X = [data(waveCorrRefno).dXY data(waveCorrRefno).maxCorr data(waveCorrRefno).AlltextureDescriptors data(waveCorrRefno).AllmomentInvariants data(waveCorrRefno).Allstatxture data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
	%     Y = [data(i).dXY data(i).maxCorr data(i).AlltextureDescriptors data(i).AllmomentInvariants data(i).Allstatxture data(i).Allregiondescriptors(:,regpropInd)];  %observation sample    
		if ~isempty(data(i).xcorrResults.dXY) && ~isempty(data(waveCorrRefno).xcorrResults.dXY)
		if waveCorrRefno == wave1number
			X = [dXY1 data(waveCorrRefno).xcorrResults.AlltextureDescriptors(:,textureInd) data(waveCorrRefno).xcorrResults.AllmomentInvariants data(waveCorrRefno).xcorrResults.Allstatxture data(waveCorrRefno).xcorrResults.Allregiondescriptors(:,regpropInd)]; %reference sample
		else
			X = [data(waveCorrRefno).xcorrResults.dXY data(waveCorrRefno).xcorrResults.AlltextureDescriptors(:,textureInd) data(waveCorrRefno).xcorrResults.AllmomentInvariants data(waveCorrRefno).xcorrResults.Allstatxture data(waveCorrRefno).xcorrResults.Allregiondescriptors(:,regpropInd)]; %reference sample
	%         X = [data(waveCorrRefno).dXY data(waveCorrRefno).Allregiondescriptors(:,regpropInd)]; %reference sample
		end
	
		Y = [data(i).xcorrResults.dXY data(i).xcorrResults.AlltextureDescriptors(:,textureInd) data(i).xcorrResults.AllmomentInvariants data(i).xcorrResults.Allstatxture data(i).xcorrResults.Allregiondescriptors(:,regpropInd)];  %observation sample            
	%     Y = [data(i).dXY data(i).Allregiondescriptors(:,regpropInd)];  %observation sample            
		%also  data(i).maxCorr  data(i).AllSqDistFromCenter
		
			if nRows <= size(Y,2)
				DistMahal = 1e07;  %arbitrary large Mahal distance number. e.g. if there were no frame information returned from the current comparision wave 
			else
				DistMahal = mahal(Y,X); % Mahalanobis
			end
		else
			DistMahal = 1e07;  %arbitrary large Mahal distance number. e.g. if there were no frame information returned from the current comparision wave 
		end
		
		dXY2 = data(i).xcorrResults.dXY;
		AllSqDistFromCenter = [];
		for j = 1:size(dXY2)
            if ~isempty(data(waveCorrRefno).xcorrResults.dXY)
			CenterXY = data(waveCorrRefno).xcorrResults.dXY(1,:);
			SqDistFromCenter = (dXY2(j,2) - CenterXY(1,2))^2 + (dXY2(j,1) - CenterXY(1,1))^2;
			AllSqDistFromCenter = [AllSqDistFromCenter; SqDistFromCenter];
            else
            AllSqDistFromCenter = [AllSqDistFromCenter; NaN];
            end
		end
		meanEuclDist = round(nanmean(AllSqDistFromCenter));
%		allmeanEuclDist = [allmeanEuclDist; meanDist];
	
	%         figure; hist(DistMahal,20)
	%         figure; plot(1:length(DistMahal),DistMahal,'-ok')
	%     figure; imagesc(DistMahal); colorbar; title([num2str(waveCorrRefno) ' - ' num2str(i)])
		meanDistMahal = nanmean(DistMahal);
	%     disp(num2str(meanDistMahal));
%		allmeanDistMahal = [allmeanDistMahal; meanDistMahal];
%	end
%	figure; plot(wave2numbers,allmeanDistMahal,'-ok'); title(['mean Mahalanobis dist from ref wave' num2str(waveCorrRefno) ' vs random']); zoom yon   %TESTING
%	figure; plot(wave2numbers,allmeanEuclDist,'-ok'); title(['Eucl dist from ref wave' num2str(waveCorrRefno) ' vs random']); zoom yon;    %TESTING
end
