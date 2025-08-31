function hF = lfp_ras(varargin)
%LFP_RAS displays a raster plot of all selected trials
%Synonym for lfp_spikeAnalysis('ras', ...), except that hF is returned as
%first value for compatibility with lfp_makepasteup.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if nargout == 0
    lfp_spikeAnalysis('ras', varargin{:});
else
    [result, spikes, ntrigs, hF] = lfp_spikeAnalysis('ras', varargin{:});
end