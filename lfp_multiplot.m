function lfp_multiplot(type, varargin)
% Displays up to 6 plots in one figure, showing the same analysis for up to
% 6 different values of some parameter, e.g. alignment event.  Remaining
% arguments are determined by <type>, which must be one of:
%   '2Dpop' - 2D Population Spike Histograms within session
%   'pop' - 1D Population Spike Histograms within session
%   '2Dsessionpop' - 2D Population Spike Histograms across sessions
%   'signif' - 1D Spike Comparison Significance Level within session
%   '2Dsignif' - 2D Population Spike Comparison Significance Level within session
%   'popcomp' - 1D Population Composition Histograms within session

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $
