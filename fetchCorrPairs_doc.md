Author: James B. Ackman  
Date: 2013-06-20 14:49:49  
Tags: manual, doc, protocols, methods  

# Detect cell pairs #
The purpose of this documentation is to demonstrate how to detect significant correlations between pairs of cells/ROIs that exhibit synchronous calcium activities.

Usage:

	function [corr_pairs,useGaussEvents,win,spkLength,numres,p_val,pvalCorrMatrix,region] = fetchCorrPairs(region,sig,numres,p_val,useGaussEvents,win)

default window is `win = 0.250` (250msec) if you don't use gaussEvents for gaussian smoothing.

Example using gaussian smoothing of the signal:

	[corr_pairs,useGaussEvents,win,spkLength,numres,p_val,pvalCorrMatrix,region] = fetchCorrPairs(region,1,1000,0.01,'true');


# Fetch Corr Pair Data #

This batch script will output a data table for downstream plotting and analysis. 

	data = batchFetchCorrpairs(filelist)