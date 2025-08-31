function hF = lfp_his(varargin)
%LFP_HIS displays a peri-event time histogram of all selected trials
%Synonym for lfp_spikeAnalysis('his', ...), except that hF is returned as
%first value for compatibility with lfp_makepasteup.


%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if nargout == 0
    lfp_spikeAnalysis('his', varargin{:});
else
    [result, spikes, ntrigs, hF] = lfp_spikeAnalysis('his', varargin{:});
end