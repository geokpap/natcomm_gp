function lfp_selectByRule(rule, varargin)
%LFP_SELECTBYRULE selects trials according to an arbitrary rule.
% lfp_selectByRule(rule)
% lfp_selectByRule(rule, rulename)

%INPUTS
% <rulename> is optional and purely cosmetic: if it is not empty, it
%   provides a label for the waitbar and a value for lfp_SelectionRule.  If
%   not given, no label is written on the waitbar.   <rulename> is
%   explicitly prohibited from being the same as any option name, and if
%   the second argument is an option name, then it is interpreted as
%   invoking that option.
% <rule> is either a matrix of parameter numbers and values for submission
%   to lfp_trialHasParams, or a string containing an arbitrary Matlab
%   expression that evaluates to true for the trials to be included. When
%   using string values in rules{:,2}, the string is interpreted by
%   replacing all occurrences of:
%       'HasParams(' with 'lfp_trialHasParams(trial,'
%       'HasPerm(' with 'lfp_trialHasPerm(trial,'
%       'HasEvent(' with 'lfp_trialHasEvent(trial,'
%       'HasTiming(' with 'lfp_trialHasTiming(trial,'
%   and then submitting it to "eval" in a context where the variable
%   <trial> is iterating through all trial numbers.  lfp_lib global
%   variables  and variables assigned by lfp_getEvtIDs are also accessible.
%   The result of the eval is then assigned as the value of
%   lfp_SelectedTrials(trial), which controls which trials are included in
%   subsequent operations. If <rule> is '' then all trials are de-selected
%   and lfp_SelectionRule becomes 'false'.
%NOTES
% <rule> can also contain lfp_TrialParams, e.g.:
% lfp_selectByRule('(lfp_TrialParams{trial}(17) >= 8) && ...
%       (lfp_TrialParams{trial}(17) ~= 100)');
% This selects all trials whose parameter #17 is greater than or equal
% to 8 and is not equal to 100.
%
% <rule> can also contain arbitrary Matlab expressions, e.g.:
% lfp_selectByRule('ismember(trial, 1:10:size(lfp_TrialIndex,1))')
% selects every 10th trial.  Note that you can only access the values of
% lfp_lib global variables, variables that are assigned values in your
% lfp_getEvtIDs file, and and the special variable <trial>.
%
% SIDE EFFECTS:
%   Sets lfp_SelectionRule
%   Logs new value of lfp_SelectionRule
%
% OPTIONS:
%   'and' - restricts selection to those trials that are already selected.
%       Has no effect if <rule> is 'true', 'false', or ''.
%   'or' - adds new selection to those trials that are already selected.
%       Has no effect if <rule> is 'true', 'false', or ''.

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

if nargin < 2 || ismember(varargin{1}, {'and', 'or'})
    rulename = '';
else 
    rulename = varargin{1};
    varargin(1) = [];
end

lfp_declareGlobals;

andflag = false;
orflag = false;

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'and'
                andflag = true;
            case 'or'
                orflag = true;
            otherwise
                error('lfp_selectByRule:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_selectByRule:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if andflag && orflag
    error('lfp_selectByRule:optconflict', ...
        '''and'' and ''or'' are mutually exclusive.');
end

if isempty(rule) || isequal(rule, 'false')
    lfp_SelectionRule = 'false';
    lfp_SelectedTrials(:) = false; %#ok<*NASGU>
elseif isequal(rule, 'true')
    lfp_SelectionRule = 'true';
    lfp_SelectedTrials(:) = true;
else
    switch class(rule)
        case 'char'
            if isempty(rulename)
                selectionrule = rule;
            else
                selectionrule = rulename;
            end
            rulestring = dg_substitute(rule, ...
                'HasParams(', 'lfp_trialHasParams(trial,' );
            rulestring = dg_substitute(rulestring, ...
                'HasPerm(', 'lfp_trialHasPerm(trial,' );
            rulestring = dg_substitute(rulestring, ...
                'HasEvent(', 'lfp_trialHasEvent(trial,' );
            rulestring = dg_substitute(rulestring, ...
                'HasTiming(', 'lfp_trialHasTiming(trial,' );
        case 'double'
            selectionrule = [ 'HasParams(' mat2str(rule) ')' ];
            rulestring = 'lfp_trialHasParams(trial, rule);';
        otherwise
            error('lfp_selectByRule:badrule', ...
                'Rule should be a string or an array of numbers');
    end
    if orflag
        lfp_SelectionRule = [selectionrule ' OR ' lfp_SelectionRule];
    elseif andflag
        lfp_SelectionRule = [selectionrule ' AND ' lfp_SelectionRule];
    else
        lfp_SelectionRule = selectionrule;
    end
    if ~lfp_NoWaitbar
        if isempty(rulename)
            waitBarName = '';
        else
            waitBarName = ['selecting trials for ' rulename];
        end
    end
    if ~lfp_NoWaitbar
        hWaitBar = waitbar(0, waitBarName, ...
            'Name', 'Selecting Trials');
    end
    lfp_getEvtIDs;
    for trial = 1:length(lfp_SelectedTrials) %#ok<*NODEF>
        rulevalue = eval(rulestring);
        if andflag
            lfp_SelectedTrials(trial) = lfp_SelectedTrials(trial) ...
                && rulevalue; %#ok<*AGROW>
        elseif orflag
            lfp_SelectedTrials(trial) = lfp_SelectedTrials(trial) ...
                || rulevalue;
        else
            lfp_SelectedTrials(trial) = rulevalue;
        end
        if ~lfp_NoWaitbar
            waitbar(trial/length(lfp_SelectedTrials), hWaitBar);
        end
    end
end
lfp_log(sprintf('Applied selection: %s', lfp_SelectionRule));
if ~lfp_NoWaitbar && exist('hWaitBar', 'var')
    close(hWaitBar);
end

