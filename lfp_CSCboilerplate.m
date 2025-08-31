function [trials, filenums, window, arglist, ...
    verboseflag, getSamplesOpts] = lfp_CSCboilerplate(arglist)
% Contains generic CSC analysis boilerplate code which is substantially
% repeated at the start of every CSC analysis function.  Well, not any
% more! This function should be called immediately on entry to a calling
% function, like this:
%   function myvalue = myfuncname(varargin)
%       [trials, filenums, window, sampledata, timepts, arglist, ...
%           verboseflag, getSamplesOpts] = lfp_CSCboilerplate(varargin);
%       [ sampledata, timepts, evtmatrix, evtidx, badtrials, trials, ...
%           filenums ] = lfp_getSamples( trials, filenums, window, ...
%           getSamplesOpts{:} );
%       trials = setdiff(trials, badtrials);
%
% When calling <myfuncname>, ALL the arguments are optional, and are given
% the empty value [] if they are missing.
%INPUTS
% arglist: the value of <varargin> immediately on entry to the calling
%   function.
%OUTPUTS
% trials: the verbatim list of trials that actually went into <sampledata>.
%   <trials> is equal to <setdiff(trials, badtrials)>, where <trials> and
%   <badtrials> are the values returned by lfp_getSamples.
% filenums: masked using lfp_SelectedFiles, [] expands to all filenums.
% window: [] expands to value of lfp_XLimAll, which may itself be [].
% verboseflag: true if 'verbose' option was specified in the raw value of
%   <arglist>.  (It is an lfp_getSamples option and is thus stripped out
%   from <arglist>.)
% sampledata: the standardized sample data returned by lfp_getSamples.
% timepts: the standardized time points returned by lfp_getSamples.
% arglist: what remains of arglist after removing all the values in
%   getSamplesOpts.
% getSamplesOpts: all the optional arguments that must be passed in to
%   lfp_getSamples.
%OPTIONS
% There is no checking for meaningless options; that is the responsibility
% of the caller, which should process all the arguments returned in
% <arglist>.  However, there is special processing for lfp_getSamples
% options, which are all moved from <arglist> to <getSamplesOpts>.  In
% addition, 'verbose' sets <verboseflag> (which is false by default), and
% 'notrunc' produces special behavior (see comments above re:
% <sampledata>).
%NOTES
%	See the "switch arglist{argnum}" statement for a list of all the option
% names that are reserved for use by lfp_getSamples, which is the union of
% all the cases in that switch statement. As of 29-Sep-2014, there is still
% no implementation here for the following options:
%       'trigfunc' (not yet implemented in lfp_getSamples)
%       'evtavg' (need to consider design and efficiency issues; such notes
%           as exist are in "lfp_lib Maintenance Volume 8.doc"; I forget if
%           there are any notes about making this uniform across CSC and
%           spike analyses, but the code should probably get pulled out of
%           lfp_getSamples)
%       'evtavg2' (ditto)

%$Rev:  $
%$Date:  $
%$Author: dgibson $

if length(arglist) <3
    window = [];
else
    window = arglist{3};
    arglist(3) = [];
end
if length(arglist) < 2
    filenums = [];
else
    filenums = arglist{2};
    arglist(2) = [];
end
if length(arglist) < 1
    trials = [];
else
    trials = arglist{1};
    arglist(1) = [];
end

% <arglist> now contains what <varargin> normally would in a function
% declared like "myvalue = myfuncname(trials, filenums, window, varargin)".
getSamplesOptidx = [];
argnum = 0;
verboseflag = false;
% Mark all lfp_getSamples options and handle any that might be used in the
% body code of the calling function:
while true
    argnum = argnum + 1;
    if argnum > length(arglist)
        break
    end
    if ~ischar(arglist{argnum})
        continue
    end
    switch arglist{argnum}
        case {
                'evtavg'
                'evtavg2'
                'evtbounds'
                }
            % opts with one arg.
            getSamplesOptidx(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            getSamplesOptidx(end+1) = argnum;
        case 'notrunc'
            getSamplesOptidx(end+1) = argnum;
        case {
                'logical'
                'multitrig'
                'norefOK'
                'rmdc'
                'rmtrend'
                'rmEP'
                'shuffle'
                'trigfunc'
                }
            % opts without arg.
            getSamplesOptidx(end+1) = argnum;
        case 'verbose'
            getSamplesOptidx(end+1) = argnum;
            verboseflag = true;
    end
end
% Move stuff as specified from <arglist> to <getSamplesOpts>:
getSamplesOpts = arglist(getSamplesOptidx);
arglist(getSamplesOptidx) = [];


