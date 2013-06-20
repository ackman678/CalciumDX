%calciumdxeventDetectArtifacts
%moving average filter with window size 2 (movement artifacts result in zaxis changes involving 100% population dF change possible spanning 2 frames)
x=nt;
% x=region.traces;
ntFilt2=zeros(size(nt));
windowSize=2;
for i=1:size(x,1)
     y=filter(ones(1,windowSize)/windowSize,1,x(i,:)'); %same as rolling window mean filter of length 'windowSize'
    ntFilt2(i,:)=y';
end
% figure; imagesc(ntFilt2)
x=mean(ntFilt2,1);
[stn, decpt] = calciumdxdetartifacts(x,ntFilt2);
% [stn, decpt] = calciumdxdetartifacts_positive(x,ntFilt2);

if ~isempty(stn)
    region.artifactFrames=[];
for i=1:length(stn)
    region.artifactFrames=[region.artifactFrames; stn(i) decpt(i)];
end
end

artifactFrameInds = zeros(size(region.traces,2),1);
for i = 1:size(region.artifactFrames,1)
artifactFrameInds(region.artifactFrames(i,1):region.artifactFrames(i,2)) = 1;
end
region.artifactFrameInds = artifactFrameInds;