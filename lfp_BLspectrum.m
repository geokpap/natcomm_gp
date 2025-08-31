function BL = lfp_BLspectrum(trials, channel, evtbounds, offsets, ...
    window, varargin)
%BL = lfp_BLspectrum(trials, channel, evtbounds, offsets, window)
%
%  <BL> is a structure containing:
%   BL.sum: a column vector containing the sum of all spectral
%       observations; the mean spectrum computed as BL.sum/BL.N is exactly
%       as computed by lfp_mtspectrum(..., 'avg').
%  	BL.sumsqrs: the unbiased estimate of the variance of the spectrum is
%       computed as (BL.sumsqrs - (BL.sum).^2/BL.N)/(BL.N-1).  The standard
%       error of the mean is thus 
%           sqrt((BL.sumsqrs - (BL.sum).^2/BL.N)/((BL.N-1)*BL.N))
%   BL.N: the number of observations (windows times tapers) that produced
%       the sums.
%   BL.f: the frequency represented by each point in BL.sum and BL.sumsqrs.
%   BL.paramlist: a cell array whose first element is the date and time of
%   the beginning of computation (as a serial date number), followed by the
%   entire argument list as it appears after argument processing, followed
%   by the values of variables that affect the computation (see code for
%   details; search for string "paramlist =").
%  To compute BL over many different session fragments without having to
%  fit them all in memory at once, simply sum each of the first three
%  fields of the BL structure returned for each fragment, and then apply
%  the same formulas given above to the sums.  (This formula for variance,
%  which is quoted at
%  http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap07/means-7.html, is
%  mathematically equivalent to the sum-of-diferences-of-squares formula
%  quoted in the Matlab help for std, which can be verified by a page or
%  two of rearranging orders of summation and collecting squared terms
%  and cross-terms.)
%
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
%  the trial is used to specify the start time. The second cell works similarly to
%  specify end time.  
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
%       more than two points with either of those values will be excluded.
%  'nodataOK' - similarly to the 'norefOK' in other functions, this
%       simply skips any trials for which the specified range of
%       <evtbounds> and <offsets> does not exist; that could be because
%       the specified range would extend beyond the containing 
%       recorded segment, or because the start event and/or end event is
%       missing in the trial.  If this option is not specified, then an
%       error is raised under this condition.
%  'k', k; 'nw', nw; 'pad', N; 'rmdc', 'session', sessionname - all of
%       these function exactly as for lfp_mtspectrum.

%$Rev: 327 $
%$Date: 2014-05-29 17:22:03 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
K = 1;
maxdatasize = 2^22;
noclip = [];
nodataOKflag = false;
NW = 1.8;
padfactor = 0;
rmdcflag = false;
session = '';
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
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
        otherwise
            error('lfp_BLspectrum:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
argnum = argnum + 1;
end

[windows, winsamples] = lfp_collectDataWindows(channel, trials, ...
    evtbounds, offsets, window, noclip, nodataOKflag);
paramlist = {now, trials, channel, evtbounds, offsets, ...
    window, varargin, NW, K, padfactor, rmdcflag, ...
    lfp_SamplePeriod, lfp_FreqLim};

% compute the sums, in batches if necessary
numwin = size(windows, 1);
if numwin == 0
    error('lfp_BLspectrum:numwin', ...
        'There are no data windows.');
end
if numwin < 10
    warning('lfp_BLspectrum:numwin2', ...
        'There are only %d data windows.', numwin);
end
startwin = 1;
batchnum = 1;
batchsize = fix(maxdatasize/(winsamples*K));
numwinleft = numwin;
windone = false(size(windows,1), 1);    % to test batching
while numwinleft > 0
    thisbatchsize = min(numwinleft, batchsize);
    sampidx = zeros(thisbatchsize, winsamples);
    for winnum = 1:thisbatchsize
        abswinnum = startwin + winnum - 1;
        sampidx(winnum,:) = windows(abswinnum, 1) : windows(abswinnum, 2);
        windone(abswinnum) = true;
    end
    data = lfp_Samples{channel}(sampidx)';
    if rmdcflag
        data = data - repmat(mean(data,1),size(data,1),1);
    end
    lfp_log(sprintf( ...
        'lfp_BLspectrum NW=%d, K=%d, pad=%d, numwinleft=%d', ...
        NW, K, padfactor, numwinleft));
    [S, BL.f, S2] = BLmtspectrumc(data, [NW K], ...
        padfactor, 1/lfp_SamplePeriod, lfp_FreqLim );
    if startwin == 1
        BL.sum = sum(S, 2);
        BL.sumsqrs = sum(S2, 2);
        BL.N = thisbatchsize;
    else
        BL.sum = BL.sum + sum(S, 2);
        BL.sumsqrs = BL.sumsqrs + sum(S2, 2);
        BL.N = BL.N + thisbatchsize;
    end
    batchnum = batchnum + 1;
    numwinleft = numwinleft - thisbatchsize;
    startwin = startwin + thisbatchsize;
end
if any(~windone)
    error('oops');
end
BL.paramlist = paramlist;
% Adjust N to account for the number of tapers summed by BLmtspectrumc:
BL.N = K * BL.N;

function [S,f,S2]=BLmtspectrumc(data,tapers,pad,Fs,fpass)
% This is identical to the chronux chronux_mtspectrumc, except:
%  1. There are no <err> or <trialiave> arguments.
%  2. S is the sum across tapers of the power, not the mean.
%  3. Instead of Serr, it returns S2, which is the sum across tapers of the
%  square of the power.

if nargin < 1; error('Need data'); end;
if nargin < 2; tapers=[3 5]; end;
if nargin < 3;pad=0;end;
if nargin < 4; Fs=1; end;
if nargin < 5; fpass=[0 Fs/2]; end;
if nargin < 6; err=0; end;
if nargin < 7; trialave=0; end;
if isempty(tapers); tapers=[3 5]; end;
if isempty(pad);pad=0;end;
if isempty(Fs); Fs=1; end;
if isempty(fpass); fpass=[0 Fs/2]; end;
if isempty(err); err=0; end;
if isempty(trialave); trialave=0;end;

[N,C]=size(data);
nfft=2^(nextpow2(N)+pad);
[f,findx]=chronux_getfgrid(Fs,nfft,fpass); 
tapers=chronux_dpsschk(tapers,N)/sqrt(Fs); % check tapers
J=chronux_mtfftc(data,tapers,nfft);
J=J(findx,:,:);
S=squeeze(sum(conj(J).*J,2));
S2=squeeze(sum((conj(J).*J).^2,2));

