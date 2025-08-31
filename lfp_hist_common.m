function histogram2D = lfp_hist_common(mode, trials, filenums, varargin)
% function lfp_hist_common(mode, trials, filenums) 
% mode:
%   1: lfp_hist_joey
%   2: lfp_hist_joey2
%   3: lfp_hist_joey3
%   4: lfp_hist_theresa
%   5: lfp_perf_t
%      lfp_hist_common(..., 'perf', n) sets number of trials to group by in
%      performance graphs 

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

figflag = false;
if nargout == 0
    figflag = true;
else
    histogram2D = [];
end

% need to add a script here to run before the rest of this fn executes, in order
% to set the BDoutput component values depending on the BDFormatNumber

if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

session = '';
perfblocksize = 10; %default the performance block size to 10 here
argnum = 1;
calflag = false;
eyeflag = false;
joyflag = false;
avgflag = false;
rectify = false;
histflag = false;
perfflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'cal'
            calflag = true;
        case 'eye'
            eyeflag = true;
        case 'joy'
            joyflag = true;    
        case 'avg'
            avgflag = true;
            trialinfo = [];
        case 'rect'
            rectify = true;
        case 'hist'
            histflag = true;
            trialinfo = [];
        case 'session'
            argnum = argnum + 1;
            session = varargin{argnum};
        case 'perf'
            perfflag = true;
            argnum = argnum + 1;
            perfblocksize = varargin{argnum};
        otherwise
            error('lfp_hist_common:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if strcmp(class(trials), 'char')
    trials = lfp_parseTrialStr(trials, session);
end
% Apply lfp_SelectedTrials as a mask to trials:
trials = trials(find(lfp_SelectedTrials(trials)));

if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_hist_common:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
            num2str(filenums(find( ...
            ~ismember(filenums, lfp_ActiveFilenums) ))) ]);
end
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_hist_common:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_hist_common:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials(find( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 ))) ]);
    % auto-indent does not work here, resetting with comment
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_hist_common:badTrials2', '<trials> must be an integer vector.');
end

switch mode
    case 1
        histogram2D = lfp_hist_j(trials, filenums, avgflag, histflag, figflag, calflag);    
    case 2
        lfp_hist_j2(trials, avgflag, histflag, calflag);    
    case 3
        lfp_hist_j3(trials, filenums, avgflag, histflag, calflag);    
    case 4
        lfp_hist_t(trials, filenums, avgflag, histflag, figflag, calflag);
    case 5
        lfp_perf_t2(trials, perfflag, perfblocksize);
end


