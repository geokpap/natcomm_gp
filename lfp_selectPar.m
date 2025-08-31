function lfp_selectPar(params_values)
%LFP_SELECTPAR selects trials based on the values of trial parameters in
% lfp_TrialParams.
%lfp_selectPar(params_values)
% Selects each trial that has ALL of the specified values for the specified
% params.  See lfp_trialHasParams for description of <params_values>.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

for trial = 1:size(lfp_TrialIndex, 1)
    lfp_SelectedTrials(trial) = lfp_trialHasParams(trial, params_values);
end