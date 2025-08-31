function SNR = dg_SE_snr(waveforms, alignpt)
% Single Electrode signal-to-noise ratios.
%INPUTS
% waveforms: one row per trigger, one column per sample.
% alignpt: the sample on which all waveforms are aligned.
%OUTPUTS
% SNR: Signal to noise ratio of each waveform.
%NOTES
% Computes standard deviation of pre-trigger samples in <waveforms> across
% all waveforms, uses that as the "noise" level, and returns ratio of
% peak-to-peak amplitude on each row of <waveforms> to noise SD.
% "Pre-trigger" is defined here as the first 5/8 of the samples up to
% <alignpt> (which is typically 8).

%$Rev: 277 $
%$Date: 2021-08-20 16:25:43 -0400 (Fri, 20 Aug 2021) $
%$Author: dgibson $

PTrange = 1 : round((5/8) * alignpt);
SD = std(reshape(waveforms(:, PTrange), [], 1));
pkval = max(waveforms, [], 2);
valval = min(waveforms, [], 2);
SNR = (pkval - valval) / SD;

