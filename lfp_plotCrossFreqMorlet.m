function [hF, hA, imagedata] = lfp_plotCrossFreqMorlet(XFMdata, varargin)
%hF = lfp_plotMorletgram(XFMdata)
%INPUTS
% XFMdata: as returned by lfp_crossFreqMorlet.
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
% Example: plot a crossFreqMorlet covariance matrix of trial 100:
%   lfp_plotCrossFreqMorlet(lfp_crossFreqMorlet(100), 'linlog');
% The value displayed for "win" is actually <offsets>.

%$Rev: 365 $
%$Date: 2015-09-29 19:57:36 -0400 (Tue, 29 Sep 2015) $
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

plotdata.align = XFMdata.align;
plotdata.figtype = 'Morlet Cross-Frequency Covariance';
plotdata.win = XFMdata.offsets;
plotdata.mode = 'gram';
plotdata.ntrigs = XFMdata.ntrigs;
plotdata.filenames = XFMdata.filenames;
switch lfp_TrialStyle
    case 'trialnums'
        plotdata.trials = dg_canonicalSeries(XFMdata.trials);
    case 'rule'
        plotdata.trials = XFMdata.trialslabel;
    otherwise
        plotdata.trials = dg_canonicalSeries(XFMdata.trials);
end
if length(XFMdata.filenames) > 1
    xfname = XFMdata.filenames{2};
else
    xfname = XFMdata.filenames{1};
end
yfname = XFMdata.filenames{1};
if linlogflag
    plotdata.ylab = sprintf('%s log10(Freq, Hz)', yfname);
    plotdata.yval = log10(XFMdata.f);
    plotdata.xlab = sprintf('%s log10(Freq, Hz)', xfname);
    plotdata.xval = log10(XFMdata.f);
else
    plotdata.ylab = sprintf('%s Freq, Hz', yfname);
    plotdata.yval = XFMdata.f;
    plotdata.xlab = sprintf('%s Freq, Hz', xfname);
    plotdata.xval = XFMdata.f;
end
%     xlabel(hA, [lfp_FileNames{2)} ' Freq, Hz']);
%     ylabel(hA, [lfp_FileNames{1)} ' Freq, Hz']);
plotdata.cbarlabel = 'cross-covariance';

    plotdata.append = sprintf( ...
        ' evtbounds %s', dg_thing2str(XFMdata.evtbounds) );

plotdata.extraline = '';
if ~isempty(varargin)
    if ~isempty(plotdata.extraline)
        delim = ' ';
    else
        delim = '';
    end
    plotdata.extraline = [plotdata.extraline delim dg_thing2str(varargin)];
end

imagedata = XFMdata.matrix;

[hF, hA] = lfp_plot(plotdata, imagedata);

