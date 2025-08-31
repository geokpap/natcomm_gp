function lfp_selectPerm(params, values)
%LFP_SELECTPERM selects trials based on any permutation of values of trial 
% parameters in lfp_TrialParams.
%lfp_selectPerm(params, values)

% see lfp_trialHasPerm for details.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

if nargin ~= 2
    error('lfp_selectPerm:nargs', 'You must supply exactly 2 arguments.');
elseif length(params) ~= length(values)
    error('lfp_selectPerm:arglengths', 'The arguments must be the same length.');
end

lfp_declareGlobals;

lfp_SelectedTrials(:) = false;
for trial = 1:size(lfp_TrialIndex, 1)
    lfp_SelectedTrials(trial) = lfp_trialHasPerm(trial, params, values);
end