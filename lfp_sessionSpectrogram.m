function [spectra, plotdata] = lfp_sessionSpectrogram(filenums, ...
    window, varargin)
%[spectra, plotdata] = lfp_sessionSpectrogram(filenums, window)
% Calls lfp_mtspectrum for each enabled trial and returns the entire set of
% results.
%INPUTS
% filenums: as for lfp_mtspectrum.
% window: relative to lfp_AlignmentRef
%OUTPUTS
% plotdata: struct from lfp_getPlotdata; fields include:
%   xvals - frequency values for the rows of <spectra>.
% spectra: spectra for each selected trial, in freqs X trials format.
%OPTIONS
% All options are passed through to lfp_mtspectrum unmodified, except for
% the following which only affect lfp_sessionSpectrogram:
%
%NOTES
% Does not plot anything, just does the computation. Plot with e.g.
%   lfp_plotSessionSpectrogram(spectra, plotdata)

%$Rev: 372 $
%$Date: 2016-01-12 18:10:51 -0500 (Tue, 12 Jan 2016) $
%$Author: dgibson $

global lfp_SamplePeriod

commontimeopts = {'recseg'};

opts2delete = [];
argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'multiwin'
                argnum = argnum + 1;
            case 'arrows'
                argnum = argnum + 1;
            case 'lines'
                argnum = argnum + 1;
            case 'avg'
            case 'EP'
            case 'err'
                argnum = argnum + 1;
            case 'evtbounds'
                argnum = argnum + 1;
                commontimeopts = [commontimeopts ...
                    {'evtbounds' varargin{argnum}}]; %#ok<*AGROW>
            case 'k'
                argnum = argnum + 1;
            case 'log'
            case 'multitrig'
                commontimeopts = [commontimeopts {'multitrig'}];
            case 'noplot'
                axinfo.plotflag = false;
            case 'notrunc'
                commontimeopts = [commontimeopts {'notrunc' window}];
            case 'norefOK'
                commontimeopts = [commontimeopts {'norefOK'}];
            case 'norm'
            case 'nw'
                argnum = argnum + 1;
            case 'coh'
            case 'ovr'
            case 'p'
                argnum = argnum + 1;
            case 'pad'
                argnum = argnum + 1;
            case 'phi'
            case 'print'
            case 'rotate'
            case 'rmBL'
                argnum = argnum + 1;
            case 'rmdc'
            case 'rmtrend'
            case 'rmEP'
            case 'session'
                argnum = argnum + 1;
            case 'showtrialnums'
            case 'thresh'
            case 'unwrap'
            case 'xscale'
                argnum = argnum + 1;
            case 'yscale'
                argnum = argnum + 1;
            otherwise
                error('lfp_sessionSpectrogram:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_sessionSpectrogram:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
varargin(opts2delete) = [];

trials = lfp_enabledTrials;
[interval, triginfo] = lfp_findCommonTime(trials, commontimeopts{:});
if any(triginfo(:,3) == 0)
    error('lfp_sessionSpectrogram:noref', ...
        'There was no reference event in trial(s) %s', ...
        dg_thing2str(dg_canonicalSeries(trials(triginfo(:,3) == 0))));
end
if interval(1) * lfp_SamplePeriod > window(1) + lfp_SamplePeriod/2 ...
        || interval(2) * lfp_SamplePeriod < window(2) - lfp_SamplePeriod/2
    error('lfp_sessionSpectrogram:window', ...
        'The set of trials selected cannot accommodate <window>.');
end
for trialidx = 1:length(trials)
    [~, ~, f, ~, ~, Spectra2plot] = lfp_mtspectrum(trials(trialidx), ...
        filenums, window, 'noplot', varargin{:});
    if ~exist('spectra', 'var')
        spectra = NaN(length(f), length(trials));
    end
    spectra(:, trialidx) = Spectra2plot(:,1);
end

plotdata = lfp_getPlotdata('sample', trials, filenums, window, triginfo, f);

