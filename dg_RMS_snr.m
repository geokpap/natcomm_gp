function SNR = dg_RMS_snr(waveforms, alignpt)
% Single Electrode signal-to-noise ratios.
%INPUTS
% waveforms: one row per trigger, one column per sample.
% alignpt: the sample on which all waveforms are aligned.
%OUTPUTS
% SNR: Signal to noise ratio of each waveform.
%NOTES
% This is the same as dg_SE_snr, except that the SD formula has been
% replaced with the simpler peak/(sqrt(2)*RMS).
%   It is impossible for SNR to be Inf unless the "baseline period"
% contains only zeros.

%$Rev: 282 $
%$Date: 2021-11-02 18:43:51 -0400 (Tue, 02 Nov 2021) $
%$Author: dgibson $

PTrange = 1 : round((5/8) * alignpt);
baselinevals = reshape(waveforms(:, PTrange), [], 1);
RMS = sqrt(mean(baselinevals.^2));
pkval = max(waveforms, [], 2);
SNR = pkval / (sqrt(2) * RMS);

