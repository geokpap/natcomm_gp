function [corr, f, windows, CV, hI, hCB] = lfp_crossFreqCorr(trials, ...
    channels, evtbounds, offsets, window, varargin)
%[corr, f, windows] = lfp_crossFreqCorr(trials, channels, evtbounds, ...
%    offsets, window)
% Collects data in the same fashion as lfp_BLspectrum, then submits it to
% dg_crossFreqCorr for analysis. The spectra are always computed using
% 'detrend' on each window.

%INPUTS
%  <trials>: Each enabled trial (i.e. lfp_SelectedTrials(trial) is 'true'
%  and <trial> is not a member of lfp_BadTrials) in the integer vector
%  <trials> is used to compute <BL>. <trials> can be a row vector or a
%  column vector, or it can be empty, in which case all selected trials are
%  used (equivalent to specifying <trials> as
%  lfp_enabledTrials(find(lfp_SelectedTrials))). If <trials> is a string,
%  then it is interpreted as containing Unique Trial IDs, i.e. the
%  combination of a session name and a trial number as it was originally
%  numbered in that session.  The session name component is optional (see
%  lfp_parseTrialStr and lfp_getTrialNum for further details).
%
%  <evtbounds> is a 1x2 cell array.  The first cell should contain
%  a list of event IDs, and the event in the list that occurs earliest in
%  the trial is used to specify the start time. The second cell works
%  similarly to specify end time.
%
%  <offsets> is a 1x2 numeric array.  Each element represents an offset in
%  seconds that is added to the time of the corresponding event in
%  <evtbounds> to compute the final start and end times that delimit the
%  data to use to compute <BL>.  The end time is truncated so that the
%  interval from start time to end time is an integral multiple of
%  <window>.
%
%  <window> specifies the width of the moving time window used to compute
%  multiple observations of the spectrum.  The step size for the moving
%  window is exactly one window, so successive windows are disjoint and
%  therefore strictly independent.
%
%OUTPUTS
%  <corr> is the cross-frequency self-correlation matrix.
%
%  <f> is the vector of frequencies represented by the matrix.
%
%  <windows> is the set of data windows collected by
%  lfp_collectDataWindows.
%
%  <CV> is the coefficient of variation, defined as standard deviation
%  divided by mean, of the power at each frequency.  It has one column for
%  each channel.
%
%  <hF> is handle to new figure window.
%OPTIONS
%   Note: 'rmEP' is not offered here because this function handles data in
%   chunks that correspond to windows, not to trials.
%  'noclip', noclip - if <noclip> is 'auto', then the maximum and minimum
%       values found in each of <channels> are treated as clipping values.
%       Otherwise, it should be a 2-element array containing the lower
%       and upper clipping values in any order.  Any windows containing
%       more than one point with either of those values will be excluded.
%  'nodataOK' - similarly to the 'norefOK' in other functions, this
%       simply skips any trials for which the specified range of
%       <evtbounds> and <offsets> does not exist; that could be because
%       the specified range would extend beyond the containing
%       recorded segment, or because the start event and/or end event is
%       missing in the trial.  If this option is not specified, then an
%       error is raised under this condition.
%  'k', k; 'nw', nw; 'pad', N; - all of
%       these function exactly as for lfp_mtspectrum.
%
% This program also recognizes and passes through to dg_crossFreqCorr these
% options:
%'freqlim', freqlim
%  <freqlim> is a 2-element row vector specifying the lower and upper
%  limits of frequencies to include in the correlation matrix.  These are
%  specified in the same units as f (see 'fs' option below).
%'k', numtapers
%  Forces number of multitapers to <numtapers>; default is 1.
%'noplot'
%  Supresses plotting of correlation matrix.
%'nw', nw
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8.
%'pad', N
%  N=0 does not pad at all (note that this differs from lfp_mtspectrum),
%  N=1 pads to twice the next greater power of 2, N=2 pads to four times
%  that length, etc.  The default value is N=0.

%$Rev: 365 $
%$Date: 2015-09-29 19:57:36 -0400 (Tue, 29 Sep 2015) $
%$Author: dgibson $

global lfp_FileNames lfp_SessionNames lfp_Samples lfp_TrialIndex ...
    lfp_FreqLim lfp_SelectionRule lfp_SamplePeriod

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

if isempty(lfp_FreqLim)
    freqlim = [0 0.5/lfp_SamplePeriod];
else
    freqlim = lfp_FreqLim;
end
K = 1;
noclip = [];
nodataOKflag = false;
NW = 1.8;
padfactor = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'k'
            argnum = argnum + 1;
            K = varargin{argnum};
        case 'noclip'
            argnum = argnum + 1;
            noclip = varargin{argnum};
        case 'nodataOK'
            nodataOKflag = true;
        case 'nw'
            argnum = argnum + 1;
            NW = varargin{argnum};
        case 'pad'
            argnum = argnum + 1;
            padfactor = varargin{argnum};
        otherwise
            error('lfp_BLspectrum:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

windows = lfp_collectDataWindows(channels(1), trials, ...
    evtbounds, offsets, window, noclip, nodataOKflag);
if length(channels) == 1
    channelstr = lfp_FileNames{channels};
else
    channelstr = dg_thing2str(lfp_FileNames(channels));
end
infostr = sprintf('%s %s trials %s sel %s evts %s offsets %s, noclip=%s, nodataOK=%d', ...
    lfp_SessionNames{1}, channelstr, ...
    dg_canonicalSeries(trials), ...
    dg_thing2str(lfp_SelectionRule), ...
    dg_thing2str(evtbounds), dg_thing2str(offsets), ...
    dg_thing2str(noclip), nodataOKflag );
if length(channels) == 1
    [corr, f, CV, hI, hCB] = dg_crossFreqCorr(windows, ...
        'getfrom', lfp_Samples{channels}, 'detrend', true, ...
        'fs', 1/lfp_SamplePeriod, ...
        'freqlim', freqlim, 'pad', padfactor, ...
        'nw', NW, 'k', K, 'descrip', infostr);
    hA = get(hI, 'Parent');
    xlabel(hA, [lfp_FileNames{channels} ' Freq, Hz']);
    ylabel(hA, [lfp_FileNames{channels} ' Freq, Hz']);
else
    [corr, f, CV, hI, hCB] = dg_crossFreqCorr(windows, ...
        'getfrom', lfp_Samples{channels(1)}, 'detrend', true, ...
        'getfrom2', lfp_Samples{channels(2)}, ...
        'fs', 1/lfp_SamplePeriod, ...
        'freqlim', freqlim, 'pad', padfactor, ...
        'nw', NW, 'k', K, 'descrip', infostr);
    hA = get(hI, 'Parent');
    xlabel(hA, [lfp_FileNames{channels(2)} ' Freq, Hz']);
    ylabel(hA, [lfp_FileNames{channels(1)} ' Freq, Hz']);
end

