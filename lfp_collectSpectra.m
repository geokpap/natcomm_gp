function [Spectra, f, windows] = lfp_collectSpectra(trials, channel, ...
    evtbounds, offsets, window, varargin)
%[corr, f, windows] = lfp_collectSpectra(trials, channel, evtbounds, ...
%    offsets, window)
% Collects data and performs batched spectral power calculation in the same
% fashion as lfp_crossFreqCorr, then concatenates the cached disk files and
% returns the result.

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
%  <Spectra> contains the power spectrum for each window, in [freq X
%  window] format. 
%
%  <f> is the vector of frequencies represented by the matrix.
%
%  <windows> is the set of data windows collected by
%  lfp_collectDataWindows.
%
%OPTIONS
%   Note: 'rmEP' is not offered here because this function handles data in
%   chunks that correspond to windows, not to trials.
%  'maxdatasize', maxdatasize - specifies the maximum number of values to
%       compute in a single batch.  (Default value of 2^22 ran
%       on my Dell Precision M65 without out-of-mem error in the full JVM
%       version for some intuitively chosen spot-test cases.  However, the
%       value required to avoid an out-of-memory error depends on the size
%       of the data loaded, as well as on whether JVM is running or not.)
%  'noclip', noclip - if <noclip> is 'auto', then the maximum and minimum
%       values found in <channel> are treated as clipping values.
%       Otherwise, it should be a 2-element array containing the lower
%       and upper clipping values in any order.  Any frames containing
%       more than one point with either of those values will be excluded.
%  'nodataOK' - similarly to the 'norefOK' in other functions, this
%       simply skips any trials for which the specified range of
%       <evtbounds> and <offsets> does not exist; that could be because
%       the specified range would extend beyond the containing 
%       recorded segment, or because the start event and/or end event is
%       missing in the trial.  If this option is not specified, then an
%       error is raised under this condition.
%  'k', k; 'nw', nw; 'pad', N; 'rmdc', 'session', sessionname - all of
%       these function exactly as for lfp_mtspectrum.
%
% This program also recognizes and passes through to dg_batchedSpectra
% these options: 
%'detrend', detrendflag
%  If <detrendflag> is true (default), applies Matlab 'detrend' function to
%  remove linear trend separately from each column of <data> before
%  processing.  Note that this implicitly includes the equivalent of 'rmdc'
%  as part of the detrending.  If false, there is no detrending or removal
%  of dc.  If <detrendflag> is 'rmdc', removes dc without removing trend.
%'freqlim', freqlim
%  <freqlim> is a 2-element row vector specifying the lower and upper
%  limits of frequencies to include in the correlation matrix.  These are
%  specified in the same units as f (see 'fs' option below).
%'k', numtapers
%  Forces number of multitapers to <numtapers>; default is 1 (note that
%  this differs from lfp_mtspectrum).
%'nw', nw
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8 (note that this differs from lfp_mtspectrum).
%'pad', N
%  N=0 does not pad at all (note that this differs from lfp_mtspectrum),
%  N=1 pads to twice the next greater power of 2, N=2 pads to four times
%  that length, etc.  The default value is N=0.

%$Rev: 76 $
%$Date: 2009-07-30 14:50:13 -0400 (Thu, 30 Jul 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

freqlim = lfp_FreqLim;
K = 3;
maxdatasize = 2^22;
noclip = [];
nodataOKflag = false;
NW = 2;
padfactor = 0;
rmdcflag = false;
session = '';
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'k'
            argnum = argnum + 1;
            K = varargin{argnum};
        case 'maxdatasize'
            argnum = argnum + 1;
            maxdatasize = varargin{argnum};
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
        case 'rmdc'
            rmdcflag = true;
        case 'session'
            argnum = argnum + 1;
            session = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_BLspectrum:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
argnum = argnum + 1;
end

[windows, winsamples] = lfp_collectDataWindows(channel, trials, ...
    evtbounds, offsets, window, noclip, nodataOKflag);
[Spectra, f] = dg_batchedSpectra(windows, ...
    'getfrom', lfp_Samples{channel}, 'detrend', true, ...
    'fs', 1/lfp_SamplePeriod, ...
    'freqlim', [0 0.5/lfp_SamplePeriod], 'pad', padfactor, ...
    'nw', NW, 'k', K);

