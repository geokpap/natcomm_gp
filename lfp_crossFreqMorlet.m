function XFMdata = lfp_crossFreqMorlet(trials, filenums, ...
    evtbounds, offsets, varargin)
%XFMdata = lfp_crossFreqMorlet(trials, filenums, ...
%    evtbounds, offsets)
% Collects data windows specified by <evtbounds> and <offsets>, then
% submits them to dg_crossFreqMorlet for analysis. The transforms are
% always computed using 'detrend' on each window.  To prevent windowing
% artifacts, the window in each trial is expanded by a margin equal to 1.5
% seconds times the largest scale used (see <margintime> in code), and the
% same amount of time is trimmed off of the scalogram before running
% statistics.
%INPUTS
% trials:  can be a row vector or a column vector, or it can be empty, in
%   which case all enabled trials are used (equivalent to specifying
%   <trials> as find(lfp_SelectedTrials)). If <trials> is a string, then it
%   is interpreted as containing Unique Trial IDs, i.e. the combination of
%   a session name and a trial number as it was originally numbered in that
%   session; the session name component is optional (see lfp_parseTrialStr
%   and lfp_getTrialNum for further details).  Only trials that are enabled
%   (i.e. lfp_SelectedTrials(trial) is 'true' and <trial> is not a member
%   of lfp_BadTrials) are actually included in the output.
% filenums:  specifies the CSC filenums to use.  If <filenums> is empty,
%   then all filenums are used.  Only filenums that are enabled (i.e.
%   lfp_SelectedFiles(filenum) is 'true') are actually included in the
%   output.
% evtbounds: see lfp_collectDataWindows; note window = 0.
% offsets: see lfp_collectDataWindows; note window = 0.
%OUTPUTS
% XFMdata: a struct suitable for use with lfp_plotXFM.  Contains the
%     following fields, the last batch of which contain the values of the
%     like-named local options variables:
%   f: column vector of pseudo-frequency for each row in the scalogram. If
%     'phase' option is used, a second column contains the phase values
%     for each column of <matrix>.
%   sum: column vector of the sum across trials of the pseudo-power of the
%       wavelet transform. Contains two columns if two <filenums> are
%       given (see lfp_BLspectrum for additional notes).
%   sumsqrs: column vector of the sum across trials of the squared
%       pseudo-power of the wavelet transform. Contains two columns if two
%       <filenums> are given.
%   matrix: cross-covariance matrix between filenums of the magnitudes
%     for each value of <f> (auto-covariance if only one channel is given).
%     Rows refer to <f> in filenums(1), columns refer to <f> in
%     filenums(2).
%   N: the number of trials that went into the 'sum' and 'sumsqrs' fields.
%   freqlim
%   noclip
%   nodataOK
%   phase
%   ntrigs: number of trials actually contributing to analysis.
%   trials: value of lfp_enabledTrials, so this might include trials that
%       were skipped pursuant to 'nodataOK'.
%   trialslabel: as computed when lfp_TrialStyle = 'rule'.
%   filenames
%   sessionnames
%   align
%   evtbounds
%OPTIONS
% 'freqlim', freqlim - <freqlim> is a 2-element row vector specifying the
%   lower and upper limits of frequencies to include in the correlation
%   matrix.  These are specified in the same units as f (see 'fs' option
%   below).  Default is lfp_FreqLim, or if that's empty, [5
%   1/(3*lfp_SamplePeriod)].
% 'fspacing', fspacing
%   ratio between scales on successive rows, which gets passed in to
%   dg_crossFreqMorlet.  Use a value of 10^(1/n) where <n> is an integer if
%   you want nice round frequency tick mark labels. Default is determined
%   by dg_crossFreqMorlet.
% 'nomatrix' - skip computation of the cross-covariance matrix, in which
%   case XFMdata.matrix is given the value [].
% 'noclip', value - <value> is used as the value for the <noclip> argument
%   to lfp_collectDataWindows.
% 'nodataOK' - if 'true', lfp_collectDataWindows issues a warning and keeps
%   going instead of raising an error when a trial has some fatal flaw in
%   the data.
% 'param', paramval - the Morlet wavelet parameter value, which gets
%   passed in to dg_crossFreqMorlet.  Default = 6.
% 'phase' - use the phase instead of the magnitude of the wavelet
%   transforms from filenums(2).
% 'verbose' - optional.
%
%EXAMPLES
%
% Ex. 1. Calculate and plot the Morlet-based crosscovariance between two
%   channels for an entire session, with full commentary:
% XFM = lfp_crossFreqMorlet([], [1 2], {lfp_NominalTrialStart ...
%       lfp_NominalTrialEnd}, [0 0], 'verbose');
% lfp_plotCrossFreqMorlet(XFM);
%
% Ex. 2. Calculate the Morlet-based coefficient of variation for
%   one channel using the last two seconds before event 8, and plot it with
%   a logarithmic frequency scale:
% XFM = lfp_crossFreqMorlet([], 1, {8 8}, [-2 0], 'nomatrix');
% SD = sqrt((XFM.sumsqrs - XFM.sum.^2/XFM.N)/(XFM.N-1));
% CV = SD ./ (XFM.sum/XFM.N);
% figure; plot(XFM.f, CV);
% set(gca, 'XScale', 'log');
% xlim(lfp_FreqLim);
%
% Ex. 3.  Fit a pink noise spectrum for use with lfp_morletgram.
% XFM = lfp_crossFreqMorlet([], 1, ...
%   {lfp_NominalTrialStart lfp_NominalTrialEnd}, [0 0], 'nomatrix');
% [pinkBL, a, b] = dg_fitPinkSpec(XFM, 'logfit');
% figure; plot(XFM.f, [XFM.sum/XFM.N pinkBL.sum])
% lfp_plotMorletgram(lfp_morletgram([], 1, [], 'norm', -b));

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

global lfp_AlignmentRef lfp_FileNames lfp_FreqLim lfp_Samples ...
    lfp_SamplePeriod lfp_SessionNames

if isempty(lfp_FreqLim)
    XFMdata.freqlim = [5 1/(3*lfp_SamplePeriod)];
else
    XFMdata.freqlim = lfp_FreqLim;
end
dg_crossFreqMorletopts = {};
matrixflag = true;
paramval = 6;
XFMdata.noclip = [];
XFMdata.nodataOK = false;
XFMdata.phase = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'freqlim'
            argnum = argnum + 1;
            XFMdata.freqlim = varargin{argnum};
        case 'fspacing'
            argnum = argnum + 1;
            dg_crossFreqMorletopts{end+1} = 'fspacing'; %#ok<AGROW>
            dg_crossFreqMorletopts{end+1} = varargin{argnum}; %#ok<AGROW>
        case 'noclip'
            argnum = argnum + 1;
            XFMdata.noclip = varargin{argnum};
        case 'nodataOK'
            XFMdata.nodataOK = true;
        case 'nomatrix'
            matrixflag = false;
        case 'param'
            argnum = argnum + 1;
            paramval = varargin{argnum};
        case 'phase'
            XFMdata.phase = true;
            dg_crossFreqMorletopts{end+1} = 'phase'; %#ok<AGROW>
        case 'verbose'
            dg_crossFreqMorletopts{end+1} = 'verbose'; %#ok<AGROW>
        otherwise
            error('lfp_crossFreqMorlet:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

FourierFactor = 4*pi/(paramval+sqrt(2+paramval^2));
margintime = 1.5 / (FourierFactor * XFMdata.freqlim(1));
windows = lfp_collectDataWindows(filenums(1), trials, ...
    evtbounds, offsets, 0, XFMdata.noclip, XFMdata.nodataOK, ...
    'margintime', margintime);
% Strictly speaking, the number of points in the margin can be +/- 1 point
% from the value used here for 'pts2trim', but we assume that makes only a
% negligible difference.
if matrixflag
    [XFMdata.f, XFMdata.sum, XFMdata.sumsqrs, XFMdata.N, XFMdata.matrix] ...
        = dg_crossFreqMorlet( ...
        windows, lfp_Samples(filenums), 'detrend', true, ...
        'fs', 1/lfp_SamplePeriod, 'freqlim', XFMdata.freqlim, ...
        'param', paramval, ...
        'pts2trim', round(margintime/lfp_SamplePeriod), ...
        dg_crossFreqMorletopts{:});
else
    XFMdata.matrix = [];
    [XFMdata.f, XFMdata.sum, XFMdata.sumsqrs, XFMdata.N] = ...
        dg_crossFreqMorlet( ...
        windows, lfp_Samples(filenums), 'detrend', true, ...
        'fs', 1/lfp_SamplePeriod, 'freqlim', XFMdata.freqlim, ...
        'param', paramval, ...
        'pts2trim', round(margintime/lfp_SamplePeriod), ...
        dg_crossFreqMorletopts{:});
end
XFMdata.ntrigs = size(windows,1);
XFMdata.trials = lfp_enabledTrials;
XFMdata.trialslabel = lfp_getTrialsLabel(trials, 'rule');
XFMdata.filenames = lfp_FileNames(filenums);
XFMdata.sessionnames = lfp_SessionNames;
XFMdata.align = lfp_AlignmentRef;
XFMdata.evtbounds = evtbounds;
XFMdata.offsets = offsets;

