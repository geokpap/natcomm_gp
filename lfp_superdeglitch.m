function [result, units] = lfp_superdeglitch(filenum, thresh, maxpts)
% Dumb wrapper to use dg_superdeglitch with lfp_createWave.
%NOTES
% Using this for its intended purpose, i.e.
%   lfp_createWave(@lfp_superdeglitch, filenum, thresh, maxpts)
% requires a lot of memory.  Specifically, your running Matlab image will
% increase in size by about twice the size of the CSC channel you are
% processing.  That can be problem for some files running on some machines.
% A different strategy could be used, which is to do the interpolations on
% the CSC channel in situ (i.e. to directly destructively edit
% lfp_Samples{filenum}), but that goes outside the scope of lfp_createWave
% and will require a new function to be written (not a big project, but I
% don't have time to do it right now as I write).

%$Rev: 339 $
%$Date: 2015-03-12 21:05:40 -0400 (Thu, 12 Mar 2015) $
%$Author: dgibson $

global lfp_Samples lfp_SamplesUnits

result = dg_superdeglitch(lfp_Samples{filenum}, thresh, maxpts);
units = lfp_SamplesUnits{filenum};
