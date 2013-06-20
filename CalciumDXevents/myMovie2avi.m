function myMovie2avi(A,fnm,graylevels)
%PURPOSE -- Save a motion jpeg avi file from a matlab array
%Need the a movie (A or A2) or binary array (A3) returned from for example wholeBrain_segmentation or wholeBrain_kmeans.m
%Need a filename, fnm
%Need to specify number of graylevels if not a binary array (e.g. 256 gray levels for a dF/F signal movie)
%USAGE -- myMovie2avi(A3,fnm)
%James B. Ackman 2/20/2013

if nargin < 3 || isempty(graylevels), graylevels = 8; end %default graylevels is good for a converted binary signal.  %change to 256 for 8 bit signals.

if nargin < 2 || isempty(fnm),
	fnm2 = [datestr(now,'yyyymmdd-HHMMSS') '.avi']; 
else
	fnm2 = [fnm(1:length(fnm)-4) datestr(now,'yyyymmdd-HHMMSS') '.avi']; 
end


%Transform the binary array into a matlab specific 'movie' data structure that can be written as an motion JPEG .avi to disk.
for fr=1:size(A,3)
	I=mat2gray(A(:,:,fr));
	[I2, map] = gray2ind(I, 8); %figure; imshow(I2,map)
	M(fr) = im2frame(I2,map);
end

%write the motion JPEG .avi to disk using auto-generated datestring based filename
vidObj = VideoWriter(fnm2)
open(vidObj)
for i =1:numel(M)
	writeVideo(vidObj,M(i));
end
close(vidObj)
