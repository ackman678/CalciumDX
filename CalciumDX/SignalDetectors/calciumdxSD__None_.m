function [onsets, offsets, param] = calciumdxSD__None_(fname,region)

onsets = cell(1,length(region.contours));
offsets = cell(1,length(region.contours));
param = [];