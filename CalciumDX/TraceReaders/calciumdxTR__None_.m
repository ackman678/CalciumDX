function [tr, trhalo, param] = calciumdxTR__None_(fname,region)
% Program used by calciumdx
% Does not read traces

param = [];
tr = zeros(length(region.contours),0);
trhalo = [];