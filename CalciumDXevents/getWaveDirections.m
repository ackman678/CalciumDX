function region = getWaveDirections(region)

%getWaveDirections.m
%--wave direction----------------------------------------------------------
%James Ackman, 1/14/2011
%updated to work by location, 5/11/2011
%polar plot example for wind speeds
%get angle from arctan (atan2) or cart2pol
% theta=atan2(Y,X)

% region.orientation.description = 'anteriormedial point = (Y,X)';   (Y,X) = row, cols
% region.orientation.value = [512 0];
locationMarkers = unique(region.location);
regionsAll = splitRegion(region);
for locationIndex = locationMarkers
    tmpregion = regionsAll{locationIndex}.region;
    results = getWaveDirectionsByLocation(tmpregion,locationIndex);
    region.wavedata{locationIndex}.waveprops.wavedirection=results.wavedirection;
end

function results = getWaveDirectionsByLocation(region, locationIndex)
disp(['locationIndex==' num2str(locationIndex)])
% if isfield(region,'orientation')
%     def_ans1 = region.orientation.description;
%     def_ans2 = region.orientation.value;
%     prompt = {'Are the following origin coordinates for the analysis of direction in the image correct?','origin coordinate in image [Y]:', 'origin coordinate in image [X]:'};
%     dlg_title = 'Input experimental parameters';
%     num_lines = 1;
%     def = {def_ans1,num2str(def_ans2(1)),num2str(def_ans2(2))};
%     answer = inputdlg(prompt,dlg_title,num_lines,def);
%     
%     results.wavedirection.origin_description = answer{1};
%     results.wavedirection.origin_coord = [str2double(answer{2}) str2double(answer{3})];
%     
%     rowdifference = results.wavedirection.origin_coord(1);
%     coldifference = results.wavedirection.origin_coord(2);
% else
%     rowdifference = 0;
%     coldifference = 0;
%     
%     def_ans1 = 'anteriormedial point = (Y,X); (Y,X) = row, cols';
%     def_ans2 = [rowdifference coldifference];
%     prompt = {'Provide description of orientation coordinate in image:','orientation coordinate in image [Y]:', 'orientation coordinate in image [X]:'};
%     dlg_title = 'Input experimental parameters';
%     num_lines = 1;
%     def = {def_ans1,num2str(def_ans2(1)),num2str(def_ans2(2))};
%     answer = inputdlg(prompt,dlg_title,num_lines,def);
%     
%     results.wavedirection.origin_description = answer{1};
%     results.wavedirection.origin_coord = [str2double(answer{2}) str2double(answer{3})];
%     
%     rowdifference = results.wavedirection.origin_coord(1);
%     coldifference = results.wavedirection.origin_coord(2);
% end

results.wavedirection.origin_description = region.orientation.description;
results.wavedirection.origin_coord = [0 0];

rowdifference = 0;
coldifference = 0;

%The following is for if measurements will be setup to use same coordinate system as the image for both hemispheres
%In this case the origin with respect to midline will be different between the two hemispheres resulting in flipped measurements with respect to midline
if coldifference == 0;
    %     rowdifference = 512;
    rowdifference = size(region.image,1); %if left of image is oriented towards front of animal
    disp('Coordinates setup with anterior to left (180degs), posterior to right (0degs)')
end

if isfield(region.wavedata{locationIndex},'wavecentroids')
    
%     results.wavedirection.theta_radians_iterativeDir={};
%     results.wavedirection.theta_radians_iterativeRho={};
    results.wavedirection.theta_radians=[];
    results.wavedirection.rho_pixels=[];
%     wfr=1;
    for wfr=1:length(region.wavedata{locationIndex}.waveonsets)
        disp(wfr)
        wavecentroid=region.wavedata{locationIndex}.wavecentroids{wfr};
        ind=1:size(wavecentroid,1);
        goodind=ind(find(~isnan(wavecentroid(:,1))));
        
        if length(goodind) > 1
            goodwavecentroids=[abs(wavecentroid(goodind,1)-coldifference) abs(wavecentroid(goodind,2)-rowdifference)];
%             goodwavecentroids=[abs(wavecentroid(goodind,1)) abs(wavecentroid(goodind,2))];
            
            % %this works--------------------------------------------------
            % for i=2:size(goodwavecentroids,1)
            % theta=atan2(goodwavecentroids(i,2)-goodwavecentroids(1,2),goodwavecentroids(i,1)-goodwavecentroids(1,1));
            % theta*(180/pi)
            % rho1 = sqrt(abs(goodwavecentroids(i,2)-goodwavecentroids(1,2)).^2 + abs(goodwavecentroids(i,1)-goodwavecentroids(1,1)).^2)
            % [theta2, rho2]= cart2pol(goodwavecentroids(i,1)-goodwavecentroids(1,1),goodwavecentroids(i,2)-goodwavecentroids(1,2));
            % theta2*(180/pi)
            % rho2
            % end
            
            % %this vectorized code also works-----------------------------
            % theta=atan2(goodwavecentroids(2:end,2)-goodwavecentroids(1,2),goodwavecentroids(2:end,1)-goodwavecentroids(1,1));
            % theta.*(180/pi)
            
            %-----Also works, the good one---------------------------------
            %Add id. direction vectors for ea wave together to find the mean direction. Calculate wave angle from this direction vector.
            
            %****try one of these next two lines--***
            % [theta, rho]= cart2pol(sum(goodwavecentroids(2:end,1)-goodwavecentroids(1,1)),sum(goodwavecentroids(2:end,2)-goodwavecentroids(1,2)));  %vector sum calculated from all direction vector differences from the first goodwavecentroid waveframe (directions from waveonset)
            [theta, rho]= cart2pol(sum(goodwavecentroids(2:end,1)-goodwavecentroids(1:end-1,1)),sum(goodwavecentroids(2:end,2)-goodwavecentroids(1:end-1,2)));  %vector sum calculated from all direction vector differences (iterative direction vector changes)
            %     theta.*(180/pi)
            
%             [x,y] = pol2cart(theta,rho);
%             figure; compass(x,y)
%             title(['vector sum wave no. ' num2str(wfr)])

            %------New one using slope estimate for direction calc---------
            %{
            linearCoef = polyfit(goodwavecentroids(:,2),goodwavecentroids(:,1),1);  %linear fit to the centroid data to find slope that matches direction of wave progression
            linearFit = polyval(linearCoef,goodwavecentroids(:,2)); %fit to one of the centroid location dimensions, may need to switch script to opposite dimension or write auto estimate algorithm of which dimension best describes the wave direction progression
%             figure; imshow(A1(:,:,wfr));
            figure; imshow(region.image);
            hold on;
            plot(goodwavecentroids(:,1),goodwavecentroids(:,2),'-', linearFit,goodwavecentroids(:,2),'r-');
            plot(goodwavecentroids(1,1),goodwavecentroids(1,2),'og');

            pxdist=sqrt((goodwavecentroids(:,1)-linearFit(1)).^2 + (goodwavecentroids(:,2)-goodwavecentroids(1,2)).^2);  %find distances from the fitted line points to the first centroid
            [val1,minInd] = min(pxdist);  %find index of line point that is a minimal distance from the 1st centroid
            [val2,maxInd] = max(pxdist);  %find index of line point that is a maximal distance from the 1st centroid
            plot(linearFit(minInd),goodwavecentroids(minInd,2),'oy');
            plot(linearFit(maxInd),goodwavecentroids(maxInd,2),'oy');
            hold off
            [theta, rho]= cart2pol(linearFit(maxInd)-linearFit(minInd),goodwavecentroids(maxInd,2)-goodwavecentroids(minInd,2)); %calculate direction of line from the minIndex fitted point to the maxInd fitted point
%}
            
            results.wavedirection.theta_radians=[results.wavedirection.theta_radians theta];
            results.wavedirection.rho_pixels=[results.wavedirection.rho_pixels rho];
            
            %         [theta2, rho2]= cart2pol(goodwavecentroids(2:end,1)-goodwavecentroids(1,1),goodwavecentroids(2:end,2)-goodwavecentroids(1,2));
            % %         theta2.*(180/pi)
            %
            %         [x,y] = pol2cart(theta2,rho2);
            %         figure; compass(x,y)
            %         title(['wave no. ' num2str(wfr)])
            %
            %     [x2,y2] = pol2cart(mean(theta2),mean(rho2));
            %     figure; compass(x2,y2)
            %     title(['mean theta rho wave no. ' num2str(wfr)])
            %
            %     [theta2, rho2]= cart2pol(mean(goodwavecentroids(2:end,1)-goodwavecentroids(1,1)),mean(goodwavecentroids(2:end,2)-goodwavecentroids(1,2)));
            %     theta2.*(180/pi)
            %
            %     [x2,y2] = pol2cart(theta2,rho2);
            %     figure; compass(x2,y2)
            %     title(['mean vector wave no. ' num2str(wfr)])
        else
            results.wavedirection.theta_radians=[results.wavedirection.theta_radians NaN];
            results.wavedirection.rho_pixels=[results.wavedirection.rho_pixels NaN];
            disp('not enough detected center of masses')
            %             region.wavespeeds_umpersec{wfr}=NaN;
        end
    end
    if isfield(region.wavedata{locationIndex}.waveprops,'wavedistance_px')
        theta=results.wavedirection.theta_radians(~isnan(results.wavedirection.theta_radians));
%         theta=(region.waveprops.wavedirection.theta_radians(~isnan(region.waveprops.wavedirection.theta_radians)))+pi;
        rho=region.wavedata{locationIndex}.waveprops.wavedistance_px(~isnan(results.wavedirection.theta_radians));
        rho = rho .* region.spaceres;
        
        %should add the following lines (need to be edited) to flip coord system (for display purposes) of the plots so that 0,90,180,270 are pointing towards post,lateral,anterior,medial directions for both hemsipheres respectively. Use same lines as from batchFetchWaveProp.m and myMakeMultiWaveSummaryPlot.m
%         centr = goodwavecentroids(1,:);
%         disp(['centr = ' num2str(centr)])
%         disp(['rowdiff = ' num2str(rowdifference)])
%         if results.wavedirection.origin_coord(1) > centr(2)
%               theta = theta .* -1;
%         end
        
        [x,y] = pol2cart(theta,rho);
        figure; compass(x,y)
        title('all waves, wavedistance in um')
    end
    
    %{
theta=atan2(wavecentroid(goodind(4),2)-wavecentroid(goodind(1),2),wavecentroid(goodind(4),1)-wavecentroid(goodind(1),1))
theta*(180/pi)

[theta, rho]= cart2pol(wavecentroid(goodind(4),2)-wavecentroid(goodind(1),2),wavecentroid(goodind(4),1)-wavecentroid(goodind(1),1))

[theta, rho]= cart2pol(wavecentroid(goodind(4),2)-wavecentroid(goodind(1),2),wavecentroid(goodind(4),1)-wavecentroid(goodind(1),1))


theta=atan2(wavecentroid(goodind(1),2)-wavecentroid(goodind(3),2),wavecentroid(goodind(1),1)-wavecentroid(goodind(3),1))
theta*(180/pi)
theta=atan2(wavecentroid(goodind(1),2)-wavecentroid(goodind(2),2),wavecentroid(goodind(1),1)-wavecentroid(goodind(2),1))
theta*(180/pi)
    %}
    
    
    % if length(goodind) > 1
    %     diffind=diff(goodind)';
    %     pxdist=sqrt((abs(diff(wavecentroid(goodind,1)))).^2 + (abs(diff(wavecentroid(goodind,2)))).^2);
    %     consecutivespeeds=(pxdist*region.spaceres)./(diffind*region.timeres)
    %     wavespeed=mean(consecutivespeeds);
    %     disp(['mean wavespeed: ' num2str(wavespeed) ' um/sec'])
    %     % pxdist=pxdist(~isnan(pxdist));
    %     region.wavespeeds_umpersec{wfr}=consecutivespeeds;
    % else
    %     disp('not enough detected center of masses')
    %     region.wavespeeds_umpersec{wfr}=NaN;
    % end
else
    disp('run getWaveCentroids.m first')
end