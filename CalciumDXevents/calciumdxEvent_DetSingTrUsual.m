function [s, d] = calciumdxEvent_DetSingTrUsual(tr,nn,period)



s=[];
d=[];
[s d] = calciumdxdettrial(tr(nn,:));
