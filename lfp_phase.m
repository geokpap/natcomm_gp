function lfp_phase(varargin)
%LFP_PHASE displays a peri-event time histogram of all selected trials
%Synonym for lfp_spikeAnalysis('phase', ..., wavefilenum)
%The wave in <wavefilenum> must have exactly one positive-going zero
%crossing per cycle.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_spikeAnalysis('phase', varargin{:});