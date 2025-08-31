function plotdata = lfp_getPlotdata(datatype, trials, channelnums, win, ...
    triginfo, xvals)
%plotdata = lfp_getPlotdata
% Creates a struct containing "boilerplate" fields that are required to
% identify how an analysis was done in lfp_lib.
%INPUTS
% trials - the literal set of trials analyzed, i.e. no unselected, bad, or
%   skipped trials.
% dataype - see 'switch datatype' for valid values.
% win - the entire time window per trial, relative to lfp_AlignmentRef,
%  that was analyzed. Analyses that use a sliding time window should
%  specify the sliding window width and step size in a separate field.
% triginfo - as returned by lfp_gatherTrialSpikes or lfp_findCommonTime.
%  Contains one row per trigger (alignment reference); column 1 =
%  start time, column 2 = end time, column 3 = trigger time; times are
%  absolute and expressed in seconds for spike data and in sample
%  indicess for sample data. Analyses that have a 'norefOK' option and
%  actually do something with trials that contain no trigger may
%  include rows where column 3 contains the value 0, indicating that
%  no trigger was found in the corresponding trial. Such rows should
%  be removed for analyses that have a 'norefOK' option and simply
%  skip the trials that are missing the trigger.
% xvals - to provide markings for x axis calibration.
%OUTPUTS
% plotdata: struct containing the following fields, which are in common to
%       all analysis routines:
%   align
%   channelnames
%   channelnums
%   sessionnames
%   trials
%   trialslabel
%   triginfo
%   win
%   xvals
%NOTES
% Any additional fields that are required should be added after calling
% lfp_getPlotdata.  Examples include the following:
%   evtavg2flag
%   evtbounds
%   multitrigflag
%   options - care should be taken not to include any large arrays here;
%       see 'trigfuncArgStr' for one method of handling them.
%   rawflag
%   slidewin - standard format is [<width> <stepsize>]
%   trigfuncArgStr = dg_thing2str(trigfuncArgs)
%   trigfuncflag
%   trigfuncH

%$Rev: 313 $
%$Date: 2013-11-27 18:00:34 -0500 (Wed, 27 Nov 2013) $
%$Author: dgibson $

global lfp_AlignmentRef lfp_SessionNames lfp_FileNames lfp_SpikeNames

plotdata.align = lfp_AlignmentRef;
plotdata.channelnums = channelnums;
plotdata.sessionnames = lfp_SessionNames;
plotdata.trials = trials;
plotdata.trialslabel = lfp_getTrialsLabel(plotdata.trials, 'rule');
plotdata.triginfo = triginfo;
plotdata.win = win;
plotdata.xvals = xvals;
switch datatype
    case 'sample'
        plotdata.channelnames = lfp_FileNames;
    case 'spike'
        plotdata.channelnames = lfp_SpikeNames;
    otherwise
        error('lfp_getPlotdata:unktype', ...
            'Unknown data type: %s', datatype);
end

