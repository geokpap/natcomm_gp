function [hF, hA, imagedata] = lfp_plotJPSTH(JPSTHdata, varargin)
%hF = lfp_plotJPSTH(JPSTHdata)
%   To make Matlab notation consistent with the notation in the original
% ref (see lfp_JPSTH), t1 has to specify the row and t2 the column, so
% channel 1 is the Y axis and channel 2 is the X axis.
%INPUTS
% JPSTHdata: as returend by lfp_JPSTH.
%OUTPUTS
% hF: figure handle to newly created figure window
% hA: axes handle containing image
% imagedata: same as the result of get(hI, 'CData') where hI is a handle to
%   the image object.
%OPTIONS
% 'diag' - plots a white line on the major diagonal.
% 'gaussthresh', plevel - thresholds the image after smoothing.  The
%   central 90% range of values after smoothing is used to estimate the
%   mean and standard deviation of a Gaussian model which is then used to
%   set threshold at the <plevel> of the Gaussian model.  Positive values
%   of <plevel> denote the right-hand tail, and negative values denote the
%   left-hand tail.  The mean of the Gaussian model is set to the mean of
%   the central 90% range of values, and the standard deviation is set to
%   1.27 times the standard deviation of the central 90% range of values.
% 'smooth', nbins - applies a Gaussian blur filter to the image having a
%   "sigma" (standard deviation) of <nbins>.  The total width of the filter
%   is 2 * <nbins> + 1 bins.
% 'thresh', percent - thresholds the image after smoothing.  If <percent>
%   is less than 50, then all values higher than <percent> percentile are
%   set to the <percent> percentile.  Otherwise, all values lower are set.
% 'trialcounts' - plots JPSTHdata.trialcountmatrix instead of
%   JPSTHdata.matrix.


%$Rev: 312 $
%$Date: 2013-10-22 16:50:43 -0400 (Tue, 22 Oct 2013) $
%$Author: dgibson $

global lfp_TrialStyle

diagflag = false;
smoothwidth = 0;
tcflag = false;
percent = [];
plevel = [];

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'diag'
                diagflag = true;
            case 'gaussthresh'
                argnum = argnum + 1;
                plevel = varargin{argnum};
            case 'smooth'
                argnum = argnum + 1;
                smoothwidth = varargin{argnum};
            case 'thresh'
                argnum = argnum + 1;
                percent = varargin{argnum};
            case 'trialcounts'
                tcflag = true;
            otherwise
                error('lfp_plotJPSTH:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_plotJPSTH:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if ~isempty(plevel) && ~isempty(percent)
    error('lfp_plotJPSTH:thresh', ...
        'You cannot use ''gaussthresh'' and ''thresh'' together.');
end

plotdata.align = JPSTHdata.align;
plotdata.figtype = 'JPSTH';
plotdata.mode = 'gram';
plotdata.ntrigs = JPSTHdata.ntrigs;
switch lfp_TrialStyle
    case 'trialnums'
        plotdata.trials = dg_canonicalSeries(JPSTHdata.trials);
    case 'rule'
        plotdata.trials = JPSTHdata.trialslabel;
    otherwise
        plotdata.trials = dg_canonicalSeries(JPSTHdata.trials);
end
plotdata.win = JPSTHdata.win;
plotdata.xlab = sprintf('%s spike time, s', JPSTHdata.spikenames{2});
plotdata.xval = JPSTHdata.timepts;
plotdata.ylab = sprintf('%s spike time, s', JPSTHdata.spikenames{1});
plotdata.yval = JPSTHdata.timepts;
if tcflag
    plotdata.cbarlabel = 'trials';
else
    if JPSTHdata.countsflag
        plotdata.cbarlabel = 'spikes^2';
    else
        plotdata.cbarlabel = 'spikes^2/sec^2';
    end
end

% binwidth etc.

if JPSTHdata.rawflag
    plotdata.figtype = [plotdata.figtype ' raw'];
end
plotdata.append = sprintf('binwidth=%d', ...
    round(1000 * median(diff(JPSTHdata.timepts))));
if JPSTHdata.evtavg2flag
    plotdata.append = [plotdata.append ' evtavg2'];
end
if ~isempty(JPSTHdata.evtbounds)
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append sprintf( ...
        ' evtbounds %s', dg_thing2str(JPSTHdata.evtbounds) )];
end
if JPSTHdata.multitrigflag
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append ' multitrig'];
end
if JPSTHdata.norefOKflag
    if ~isempty(plotdata.append)
        plotdata.append = [plotdata.append ' '];
    end
    plotdata.append = [plotdata.append ' norefOK'];
end

plotdata.extraline = '';
% info requiring plotdata.extraline
if JPSTHdata.trigfuncflag
    plotdata.extraline = sprintf('trigfunc: %s(%s)', ...
        char(JPSTHdata.trigfuncH), JPSTHdata.trigfuncArgStr);
end
if ~isempty(varargin)
    if ~isempty(plotdata.extraline)
        delim = ' ';
    else
        delim = '';
    end
    plotdata.extraline = [plotdata.extraline delim dg_thing2str(varargin)];
end
    
if tcflag
    imagedata = JPSTHdata.trialcountmatrix;
else
    imagedata = JPSTHdata.matrix;
end

if smoothwidth ~= 0
    fw = 2 * smoothwidth + 1;
    gw = repmat(normpdf((-smoothwidth:smoothwidth)/smoothwidth)', 1, fw) ...
        .* repmat(normpdf((-smoothwidth:smoothwidth)/smoothwidth), fw, 1);
    imagedata = conv2(imagedata, gw, 'same');
end

if ~isempty(plevel)
    Ltail = prctile(imagedata(:), 5);
    Rtail = prctile(imagedata(:), 95);
    ctr = imagedata > Ltail & imagedata < Rtail;
    mu = mean(imagedata(ctr));
    sigma = 1.27 * std(imagedata(ctr));
    if plevel < 0
        thresh = norminv(-plevel, mu, sigma);
        imagedata(imagedata>thresh) = thresh;
    else
        thresh = mu - sigma * norminv(plevel, 0, 1);
        imagedata(imagedata<thresh) = thresh;
    end
end

if ~isempty(percent)
    thresh = prctile(imagedata(:), percent);
    if percent < 50
        imagedata(imagedata>thresh) = thresh;
    else
        imagedata(imagedata<thresh) = thresh;
    end
end

[hF, hA] = lfp_plot(plotdata, imagedata);
if diagflag
    set(hA, 'NextPlot', 'add');
    plot([JPSTHdata.timepts(1) JPSTHdata.timepts(end)], ...
        [JPSTHdata.timepts(1) JPSTHdata.timepts(end)], 'w');
end
    
