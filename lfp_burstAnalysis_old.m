function [highscore, lowscore, highratio, plotdata] = ...
    lfp_burstAnalysis(trials, filenums, win, varargin)
%OUTPUTS
% highscore, lowscore, highratio:  as returned from dg_computeCoactivation.
% plotdata:  a structure containing the following fields:
% 	timepts:  time points for the centers of the time bins constructed
%       by dg_computeCoactivation.
%   trials:  the list of enabled trialnums actually analyzed.
%   align:  the alignment reference event for the analysis.
%   win:  goes to lfp_gatherTrialData.
%OPTIONS
% Any of the options to dg_computeCoactivation can be specified here and
%   will be passed through on the call to dg_computeCoactivation.
%NOTES
% Preserves logical data type of lfp_Samples only if ALL channels are
% logical.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It seems fishy that highratio(1,:) is always equal to either
% highratio(2,:) or highratio(3,:).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$Rev:  $
%$Date:  $
%$Author:  $

lfp_declareGlobals;

argnum = 1;
binwidth = 5;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'binsize'
            argnum = argnum + 1;
            binwidth = varargin{argnum};
        case 'clevel'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'numsets'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'numshuffles'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'offset'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'onset'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'plevel'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        case 'verbose'
            argnum = argnum + 1;
            % do nothing, just pass into dg_computeCoactivation
        otherwise
            error('lfp_burstAnalysis:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if nargin < 3
    win = [];
end

if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(lfp_SelectedFiles(filenums));
logicalflag = true;
for fileidx = 1:numel(filenums)
    if ~islogical(lfp_Samples{filenums(fileidx)}) %#ok<*USENS>
        logicalflag = false;
    end
end

if nargin < 1 || isempty(trials)
    trials = 1:length(lfp_SelectedTrials);
end
if ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
% <trials> is now numeric, i.e. trialnums.
if any(trials > length(lfp_SelectedTrials))
    warning('lfp_burstAnalysis:trials', ...
        'Ignoring trials %s, which are beyond the last trial.', ...
        dg_canonicalSeries(trials(trials > length(lfp_SelectedTrials))) );
    trials(trials > length(lfp_SelectedTrials)) = [];
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

gathertrialdataopts = {};
if logicalflag
    gathertrialdataopts{end+1} = 'logical';
end
[burstdata, interval, win] = ...
    lfp_gatherTrialData(filenums, trials, win, gathertrialdataopts{:});

[highscore, lowscore, highratio, rawhist, shufhist, rawcoact, ...
    binends] = dg_computeCoactivation_old(burstdata, varargin{:}); %#ok<*ASGLU>

rawtimepts = lfp_SamplePeriod * (interval(1):interval(2));
binnedtimepts = (rawtimepts(binends - binwidth + 1) ...
    + rawtimepts(binends)) / 2;

plotdata.timepts = binnedtimepts;
plotdata.trials = trials;
plotdata.align = lfp_AlignmentRef;
plotdata.win = win;

