function [hF, hA, imagedata] = lfp_plotMorletgram(morletdata, varargin)
%hF = lfp_plotMorletgram(morletdata)
%INPUTS
% morletdata: as returned by lfp_morletgram.
%OUTPUTS
% hF: figure handle to newly created figure window
% hA: axes handle containing image
% imagedata: same as the result of get(hI, 'CData') where hI is a handle to
%   the image object.
%OPTIONS
% 'linlog' - instead of plotting frequency on a log scale (as it is
%   furnished by lfp_morletgram), convert top and bottom freqs to log(freq)
%   and plot THAT on a linear scale for the sake of the data picker.
%NOTES
%
% Example: plot a Morlet scalogram of trial 100 (aassuming there is only
% one channel loaded):
%   lfp_plotMorletgram(lfp_morletgram(100));
%
% Example: plot a Morlet scalogram averaged over selected trials of filenum
% 3, with a time window of [-3.5 5.5] (which given the low frequency limit
% of 1 Hz will get trimmed down to [-2 4]), and normalized to a pink noise
% spectrum with an exponent of 1.9, plotted in log10 units:
%   lfp_plotMorletgram(lfp_morletgram([], 3, [-3.5 5.5], 'norm', 1.9), 'linlog')
%
% Note that the Y coordinates in the displayed image are never in Hz, even
% though the standard display shows tick marks calibrated in Hz.  However,
% for the purpose of plotting additional lines on the image, you can use
% the coordinates that are shown by the data picker, regardless of whether
% the display is standard or 'linlog'.

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

global lfp_TrialStyle

argnum = 1;
linlogflag = false;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'linlog'
                linlogflag = true;
            otherwise
                error('lfp_plotMorletgram:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_plotMorletgram:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1; 
end

% <plotdata> contains the metadata required by lfp_plot to produce the
% usual labels.
plotdata.align = morletdata.align;
plotdata.figtype = 'Morlet Scalogram';
plotdata.mode = 'gram';
plotdata.ntrigs = morletdata.ntrigs;
plotdata.filenames = morletdata.filenames;
switch lfp_TrialStyle
    case 'trialnums'
        plotdata.trials = dg_canonicalSeries(morletdata.trials);
    case 'rule'
        plotdata.trials = morletdata.trialslabel;
    otherwise
        plotdata.trials = dg_canonicalSeries(morletdata.trials);
end
plotdata.win = morletdata.win;
plotdata.xlab = sprintf('Time, s');
plotdata.xval = morletdata.timepts;
if linlogflag
    plotdata.ylab = sprintf('log10(Frequency, Hz)');
    plotdata.yval = log10(morletdata.freqs);
else
    plotdata.ylab = sprintf('Frequency, Hz');
    plotdata.yval = morletdata.freqs;
end
plotdata.cbarlabel = 'Pseudo-Power, dB';

if morletdata.evtavg2flag
    plotdata.append = [plotdata.append ' evtavg2'];
end
if ~isempty(morletdata.evtbounds)
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append sprintf( ...
        ' evtbounds %s', dg_thing2str(morletdata.evtbounds) )];
end
if morletdata.multitrigflag
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append ' multitrig'];
end
if morletdata.norefOKflag
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append ' norefOK'];
end

plotdata.extraline = morletdata.options;
% info requiring plotdata.extraline
if morletdata.trigfuncflag
    plotdata.extraline = sprintf('trigfunc: %s(%s)', ...
        char(morletdata.trigfuncH), morletdata.trigfuncArgStr);
end
if ~isempty(varargin)
    if ~isempty(plotdata.extraline)
        delim = ' ';
    else
        delim = '';
    end
    plotdata.extraline = [plotdata.extraline delim dg_thing2str(varargin)];
end

imagedata = 10*log10(morletdata.P);

[hF, hA] = lfp_plot(plotdata, imagedata);
if numel(morletdata.trials) == 1
    lfp_plotEvtMarkers(hA, morletdata.evtmatrix);
else
    lfp_plotEvtMarkers(hA, [0 morletdata.align(1)]);
end

