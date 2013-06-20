%bug testing 2011-12-14
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