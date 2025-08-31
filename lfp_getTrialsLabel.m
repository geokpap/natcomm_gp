function label = lfp_getTrialsLabel(trials, style)
%label = lfp_getTrialsLabel(trials, style)
%   Returns a string suitable for use in figure titles describing the
%   selection of trials that contributed to an aggregated plot.  <trials>
%   is a list of the trial numbers.  <style> controls what style of label
%   is returned, which may be:
%   'trialnums' - a string rendition of <trials>
%   'rule' - the selection rule, plus the value of lfp_BadTrials, plus an
%   indicator of whether the list of trial nums is identical to the list of
%   enabled trials ('all of') or is a 'superset', 'subset', or 'alteration'
%   (neither a superset nor a subset)

%$Rev: 313 $
%$Date: 2013-11-27 18:00:34 -0500 (Wed, 27 Nov 2013) $
%$Author: dgibson $

global lfp_SelectionRule lfp_BadTrials

if isempty(style)
    error('lfp_getTrialsLabel:badstyle2', ...
            'Trial selection label style is empty');
end
switch style
    case 'trialnums'
        label = dg_canonicalSeries(trials);
    case 'rule'
        label = lfp_SelectionRule;
        if ~isempty(label)
            label = [ ' of ' label];
        end
        if ~isempty(lfp_BadTrials)
            label = [ label ' excluding ' dg_canonicalSeries(lfp_BadTrials) ];
        end
        allenabledtrials = lfp_enabledTrials;
        if isequal(trials, allenabledtrials)
            label = [ 'all' label ];
        else
            if all(ismember(trials, allenabledtrials))
                label = [ 'subset' label ];
            elseif all(ismember(allenabledtrials, trials))
                label = [ 'superset' label ];
            else
                label = [ 'alteration' label ];
            end
        end
    otherwise
        error('lfp_getTrialsLabel:badstyle', ...
            'Style "%s" is not recognized', style );
end