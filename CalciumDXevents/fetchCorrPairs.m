function [region] = fetchCorrPairs(region,sig,numres,p_val,useGaussEvents,win)
%fetchCorrPairs Find temporally correlated cell pairs using Monte Carlo simulations
%[region] = fetchCorrPairs(region,sig,numres,p_val,useGaussEvents,win)
%
%Examples:
%[region] = fetchCorrPairs(region,1,1000,0.01,'true');
%[region] = fetchCorrPairs(region,1,1000,0.01,'false',0.1);
%
%Options:
%sig = signal for spike occurence, 1 for a spike matrix
%numres = no. of simulations
%p_val = significance level
%useGaussEvents; Whether or not to use spike signal smoothing with a gaussian. If false or empty, then a square wave/block signal will be used surrounding the spike at +/- win.
%win = +/- spike time window in seconds for determining correlations.  Only used if useGaussEvents = false;   
%
%Output:
% Will pass the region data structure that was input back to workspace with the new corr data at region.userdata.corr
%
%Versions:
%updated 2013-07-18 13:22:09 by jba with consolidated outputs inside the region.userdata.corr data structure
%%update with additional outputs James B. Ackman 2013-01-05 23:03:38
%2011-10-20 James B. Ackman-- added new variables for choosing whether to use Gauss events or not for spk smoothing, as well as window length defining.
%edited by jba on 30/07/2007 so that this script can still be used whether or not the region.userdata.schmutzr field is defined or not
%
%See also:
% batchFetchCorrPairs, myPlotCorrGraphImage, myPlotPvalueCorrMatrix

if nargin < 5 || isempty(useGaussEvents); useGaussEvents = 'false'; end  %smooth spk with Gaussian. Not very different results than just using ~7fr spk smoothing for datasets with low background spontaneous activity frequencies.
if nargin < 6 || isempty(win); win = 0.250; end  %250msec spk corr window

%all_s = cell(1,2); %jba
all_s ={}; %jba
all_s{1} = zeros(size(region.traces)); % setup spike matrix for all cells
for i = 1:size(region.traces,1)
    all_s{1}(i,region.onsets{i})=1;
end


%{
%--------------------------------------------------------------------------
%This section is for setting up all_s if you wanna do corr on the gdp, sch
%subsets of cells. Comment or uncomment if you wanna do these corrs.
if isfield(region,'userdata') %jba
if isfield(region.userdata,'schmutzon') %jba
%if isfield(region.userdata,'schmutzr') %added by jba
all_s{2} = all_s{1}; % non-SCH cells
all_s{2}(region.userdata.schmutzr,:) = [];
if isfield(region.userdata,'schmutzon') %added by jba
all_s{3} = zeros(length(region.userdata.schmutzon),size(region.traces,2)); % SCH cells
for i = 1:length(region.userdata.schmutzon)
    all_s{3}(i,region.userdata.schmutzon{i})=1;
end
end %jba
end %jba
end %jba
%--------------------------------------------------------------------------
%}

str = {'All cells ','Non-SCH cells ','SCH cells '};
% ii=1; %testing purposes
for ii = 1:length(all_s) %jba
    fprintf('\n');
    fprintf(str{ii});
    s = all_s{ii};
    
    %     %     s = zeros(size(region.traces));
    %     s = zeros(10,size(region.traces,2));
    % %     s(1,c1) = 1;
    % %     s(2,c2) = 1;
    % %     s(1,[150 250]) = 1;
    % %     s(2,[150 250]) = 1;
    %
    %
    %     s(1:10,[250]) = 1;
    % %     s(2,[250]) = 1;
    % %     s(3,[250]) = 1;
    %     desc = '';
    
    %-------Smooth the spike matrix--------------
    %smooth with gaussian
    if strcmp(useGaussEvents,'true')
        s_thick = gauss_events(s,sig); %smooth with gaussian. For sig=1, resulting spike matrix will be [0.0001 0.0183 0.3679 1.0000 0.3679 0.0183 0.0001] surrounding the spike onset.
        tmp = zeros(1,size(s,2));
        tmp(round(numel(tmp)/2))=1;
        tmpSmooth = gauss_events(tmp,sig);  %calculate spkLength, cause gauss_events is dependent on total no. of samples (frames)
        spkLength=numel(tmpSmooth(tmpSmooth>0.01));  %below this the gaussian drops to near zero-infinity, here it would be number of frames within the 99% confidence interval for the gaussian smoothed spike shape
        win = ((spkLength-1)/2)*region.timeres;
    else
        %smooth spike with +/- N frames, where N is set by winFrames
        s_thick = s;  %or for 1fr/spk, no smoothing. For low frame rate movies
        winFrames = round(win/region.timeres);
        %     winFrames=1;
        for c = 1:size(s,1)
            for f = find(s(c,:))
                in = find(abs((1:size(s,2))-f)<=winFrames);  %e.g. for winFrames of 2, this will give a 5 frame total spk length
                s_thick(c,in) = s_thick(c,f);
            end
        end
        spkLength = winFrames*2+1;
    end
    %-------Begin correlation-------------
    corrs = s_thick*s_thick';  %corr distance metric for the observed spike matrix
    corrs_obs = corrs;
    count = zeros(size(corrs));
    corrs_cont = reshape(corrs,1,numel(corrs)); %reshape adjacency matrix into a vector
    %     [r,p] = corrcoef(s_thick')
    corrs_res = zeros(fix(p_val*numres)+1,numel(corrs));  %this line uses the pval. Less rows from smaller p_vals (probability of event in no. of simulations) or less num_reshuffles.
    %     s_res = zeros(size(s,1),size(s,2),numres);
    npairs = zeros(1,numres);
    for t = 1:numres  %repeat the simulation numres times
        for c = 1:size(s,1);
            s(c,:) = s(c,randperm(size(s,2)));  %shuffle the spike times for each cell
        end
        %         s_res(:,:,t) = s;
        s_thick = gauss_events(s,sig); %smooth resulting shuffled spike matrix with a gaussian
        corrs = s_thick*s_thick'; %corr distance metric for the reshuffled spike matrix
        
        count = count + (corrs>=corrs_obs);
        corrs_res(1,:) = reshape(corrs,1,numel(corrs)); %add the reshuffled adjacency matrix corr values to the first row of probability array.
        %         pairs = (corrs_cont>corrs_res(1,:));  %jba
        corrs_res = sort(corrs_res); %sort columns of probablity matrix in order of ascending corr value
        if mod(t,10) == 0
            fprintf('.');
        end
        pairs = (corrs_cont>corrs_res(2,:));  %jba
        pairs = reshape(pairs,size(s,1),size(s,1));
        [i j] = find(pairs==1);
        f = find(j>i);
        pairs = [i(f) j(f)];
        npairs(t) = 100*size(pairs,1)/(size(s,1)*((size(s,1)-1))/2);
    end
    pvalCorrMatrix = count/numres;
    corr_pairs{ii} = pairs;
    %     myCDF(npairs,'percent pairs')
    fprintf('\n');
    
    %     figure;
    %     for j = 1:size(s_res,3)
    %         imagesc(s_res(:,:,j));
    %         F(j) = getframe;
    %     end
    %     movie(F,10)
    
    
    %{
    pairs_tmp = [];
    for idx = 1:size(corrs_res,1)
        pairs = (corrs_cont>corrs_res(idx,:)); %find those observed corr values exceeding the reshuffled simulation values above corrs_res(2,:)
        pairs = reshape(pairs,size(s,1),size(s,1));
        %         figure; imagesc(pairs)
        [i j] = find(pairs==1);
        f = find(j>i);
        pairs = [i(f) j(f)];
        %     corr_pairs{ii} = pairs;
        pairs_tmp = [pairs_tmp 100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2)];
%         fprintf(['Percent pairs correlated, corr_res ' num2str(idx) ': ' num2str(100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2)) '\n']);
    end
    figure; plot(1:size(corrs_res,1),pairs_tmp); ylabel('Percent pairs correlated'); xlabel(['corrs\_res(i,:), max = ' num2str(size(corrs_res,1))]); xlim([0 size(corrs_res,1)]);
    title(['numres = ' num2str(numres) ', pval = ' num2str(p_val) ', ' desc])
    
        pairs = (corrs_cont>corrs_res(2,:)); %find those observed corr values exceeding the reshuffled simulation values above corrs_res(2,:)
        pairs = reshape(pairs,size(s,1),size(s,1));
        [i j] = find(pairs==1);
        f = find(j>i);
        pairs = [i(f) j(f)];
        corr_pairs{ii} = pairs;
%         fprintf(['Percent pairs correlated, corr_res ' num2str(idx) ': ' num2str(100*size(pairs,1)/(size(s,1)*(size(s,1)-1)/2)) '\n']);
    %}
    
    fprintf(['        Total number of cells: ' num2str(size(s,1)) '\n']);
    fprintf(['Percent cells in correlations: ' num2str(100*length(unique(reshape(pairs,1,numel(pairs))))/size(s,1)) '\n']);
    fprintf(['        Total number of pairs: ' num2str(size(s,1)*(size(s,1)-1)/2) '\n']);
    fprintf(['     Percent pairs correlated: ' num2str(100*size(pairs,1)/(size(s,1)*((size(s,1)-1))/2)) '\n']);
    
    
    if ~isfield(region,'userdata') || ~isfield(region.userdata,'corr')
        region.userdata.corr{1}.corr_pairs=corr_pairs;
        region.userdata.corr{1}.pvalCorrMatrix = pvalCorrMatrix;
        region.userdata.corr{1}.params.useGaussEvents=useGaussEvents;
        region.userdata.corr{1}.params.window_sec=win;  % The user defined spike window for useGaussEvents = False. %For useGaussEvents = True it is number of frames within the 99% confidence interval for the gaussian smoothed spike shape
        region.userdata.corr{1}.params.spkLength_frames=spkLength;
        region.userdata.corr{1}.params.numres=numres;
        region.userdata.corr{1}.params.p_val=p_val;
        region.userdata.corr{1}.params.date = datestr(now,'yyyymmdd-HHMMSS');
    else
        len=length(region.userdata.corr);
        region.userdata.corr{len+1}.corr_pairs=corr_pairs;
        region.userdata.corr{len+1}.pvalCorrMatrix = pvalCorrMatrix;
        region.userdata.corr{len+1}.params.useGaussEvents=useGaussEvents;
        region.userdata.corr{len+1}.params.window_sec=win;
        region.userdata.corr{len+1}.params.spkLength_frames=spkLength;
        region.userdata.corr{len+1}.params.numres=numres;
        region.userdata.corr{len+1}.params.p_val=p_val;
        region.userdata.corr{len+1}.params.date = datestr(now,'yyyymmdd-HHMMSS');
    end
    
    
end

function s_thick = gauss_events(s,sig)

s_thick = zeros(size(s));
for c = 1:size(s,1)
    for d = find(s(c,:)==1)
        s_thick(c,:) = s_thick(c,:) + exp(-(((1:size(s,2))-d)/sig).^2);
    end
    %     disp(s_thick(s_thick>0))
end
