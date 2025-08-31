function [hI, hCB, stuff] = lfp_specPopAvg(trials, filenums, moving_win, ...
    varargin)
%[hI, hCB, stuff] = popavgmtspec(trials, filenums, moving_win)
% Computes a population average spectrogram across multiple CSC channels.
%INPUTS
% As per lfp_spec.
%OUTPUTS
% As per dg_avgGram.
%OPTIONS
% All options are passed through verbatim to lfp_spec except as noted here:
% 'nodisplay' - affects visibility of the figures created by dg_avgGram in
%   addition to those created by lfp_spec.
% 'rmBL', BL - <BL> must be an array with the same number of elements as
%   <filenums>.  The corresponding element of BL is used when lfp_spec is
%   called to analyze one of the <filenums>.

%$Rev: 300 $
%$Date: 2013-05-01 18:47:38 -0400 (Wed, 01 May 2013) $
%$Author: dgibson $

global lfp_SelectedFiles

if isempty(filenums)
    filenums = find(lfp_SelectedFiles);
end

argnum = find(cellfun(@isequal, varargin, repmat({'rmBL'}, ...
    size(varargin) )));
if isempty(argnum)
    BLopts = {};
else
    BL = varargin{argnum+1};
    BLopts = {'rmBL'};
end
varargin([argnum argnum+1]) = [];
if any((cellfun(@isequal, varargin, repmat({'nodisplay'}, ...
        size(varargin) ))))
    dispopt = {'nodisplay'};
else
    dispopt = {};
end

for fileidx = 1:length(filenums)
    if ~isempty(BLopts)
        BLopts{2} = BL(fileidx);
    end
    hF(fileidx) = lfp_spec('mt', trials, filenums(fileidx), ...
        moving_win, varargin{:}, BLopts{:} );
    fprintf('Created spec fig %d\n', hF(fileidx));
end

[hI, hCB, stuff] = dg_avgGram(hF, 'sem', [1 3 9], 'semwidth', 2, ...
    'std', dispopt{:});
