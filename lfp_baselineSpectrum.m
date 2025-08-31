function blspec = lfp_baselineSpectrum(trials, filenums, ...
    evtbounds, evtoffsets, varargin)
%blspec = lfp_baselineSpectrum(trials, filenums, ...
%     evtbounds, evtoffsets)
%blspec = lfp_baselineSpectrum(..., 'err')
%blspec = lfp_baselineSpectrum(..., 'hfnorm', f)
%blspec = lfp_baselineSpectrum(..., 'k', numtapers)
%blspec = lfp_baselineSpectrum(..., 'lin')
%blspec = lfp_baselineSpectrum(..., 'norm')
%blspec = lfp_baselineSpectrum(..., 'nw', bandwidth)
%blspec = lfp_baselineSpectrum(..., 'pad', N)
%blspec = lfp_baselineSpectrum(..., 'rmdc')
%blspec = lfp_baselineSpectrum(..., 'spikenorm')

%blspec = lfp_baselineSpectrum(trials, filenums, ...
%     {startIDs stopIDs}, evtoffsets)
%   <trials> and <filenums> are as for lfp_disp.  The third argument is a
%   cell array works like the argument to the 'evtbounds' option in
%   lfp_disp: <startIDs> is a list of alternative event IDs to use as the
%   start of the time range, and <endIDs> is a list of alternative event
%   IDs to use as the end of the time range.  The start event is the
%   earliest event in the trial that is a member of <startIDs>; the end
%   event is the earliest event in the trial AFTER the start event that is
%   a member of <endIDs>.  <evtoffsets>
%Um,wait,no: this should work like lfp_NewTrialTarget.

%blspec = lfp_baselineSpectrum(..., 'err')
%   Include 

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
session = '';
trialstyle = lfp_TrialStyle;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'avg'
            avgflag = true;
            trialinfo = [];
    end
    argnum = argnum + 1;
end
