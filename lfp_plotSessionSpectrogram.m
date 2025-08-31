function [hF, hI, hCB] = lfp_plotSessionSpectrogram(spectra, plotdata)
%INPUTS
% spectra, plotdata: as returned by lfp_sessionSpectrogram.
%OUTPUTS
% hF: handle to newly created figure.
% hI, hCB: as returned by dg_showGram.

%$Rev: 409 $
%$Date: 2020-04-25 15:13:16 -0400 (Sat, 25 Apr 2020) $
%$Author: dgibson $

global lfp_SelectionRule lfp_CLimAll

xlabstr = sprintf('selected trials: %s', lfp_SelectionRule);
hF = figure;
[hI, hCB] = dg_showGram(hF, 1:size(spectra,2), plotdata.xvals, spectra, ...
    sprintf('%s %s', plotdata.sessionnames{1}, plotdata.channelnames{1}), ...
    xlabstr, 'frequency, Hz', 'Power, dB');
if ~isempty(lfp_CLimAll)
    dg_recolorGram(hCB, lfp_CLimAll, hI);
end

