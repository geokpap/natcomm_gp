% lfp_setup_georgios_ApAv_rev_only
% The brute force "make it work" event codes, including codes from
% "Eventmarkers_ApAv.txt".  Works with session types 'ApAv' and 'ApAv_em'.
%
% 29-Feb-2024: Note that FORCED and FixPtOn are used regardless of stim type, so the
% stim block is distinguished by ['choiceCueOnEM' 'forcedRedOnEM'
% 'forcedYelOnEM'] in place of ['choiceCueOn' 'forcedRedOn'
% 'forcedYelOn'].
%   Also works for parsing 'ApAp' sessions, but I think some event codes
% are different for those (e.g. [220 221 235 236]).

%$Rev: 45 $
%$Date: 2023-03-31 20:14:39 -0400 (Fri, 31 Mar 2023) $
%$Author: dgibson $

lfp_getEvtIDs;
lfp_AutoassignFilenums = false;
lfp_SetupType = 'georgios_ApAv';
lfp_TrialStyle = 'rule';
lfp_AlignStyle = 'name';
lfp_EventDefaultColor = [.5 .5 .5];
lfp_NominalTrialStart = MLTrialStart;
lfp_NominalTrialEnd = MLTrialEnd;
lfp_LogFileName = 'ken.log';
lfp_NoWaitbar = true;

lfp_EventNames(1:200) = {'monkey_logic'};
lfp_EventNames{AVE_OFF_FORCED} = 'AVE_OFF_FORCED';
lfp_EventNames{AVE_ON_FORCED} = 'AVE_ON_FORCED';
lfp_EventNames{CROSS_LEFT} = 'CROSS_LEFT';
lfp_EventNames{CROSS_RIGHT} = 'CROSS_RIGHT';
lfp_EventNames{FORCED} = 'FORCED';
lfp_EventNames{FixPtOn} = 'FixPtOn';
lfp_EventNames{ITI1} = 'ITI1';
lfp_EventNames{ITI2} = 'ITI2';
lfp_EventNames{MLTrialEnd} = 'MLTrialEnd';
lfp_EventNames{MLTrialStart} = 'MLTrialStart';
lfp_EventNames{RWD_OFF_FORCED} = 'RWD_OFF_FORCED';
lfp_EventNames{RWD_ON_FORCED} = 'RWD_ON_FORCED';
lfp_EventNames{TargOffNoResp} = 'TargOffNoResp';
lfp_EventNames{airpuffOffAp} = 'airpuffOffAp';
lfp_EventNames{airpuffOffNoResp} = 'airpuffOffNoResp';
lfp_EventNames{airpuffOnAp} = 'airpuffOnAp';
lfp_EventNames{airpuffOnNoResp} = 'airpuffOnNoResp';
lfp_EventNames{badTargAcq} = 'badTargAcq';
lfp_EventNames{choiceCueOn} = 'choiceCueOn';
lfp_EventNames{choiceCueOnEM} = 'choiceCueOnEM';
lfp_EventNames{crossAcq} = 'crossAcq';
lfp_EventNames{crossOff} = 'crossOff';
lfp_EventNames{event_zero} = 'event_zero';
lfp_EventNames{fixbreak} = 'fixbreak';
lfp_EventNames{fixptOffFixbreak} = 'fixptOffFixbreak';
lfp_EventNames{forcedRedOff} = 'forcedRedOff';
lfp_EventNames{forcedRedOffEM} = 'forcedRedOffEM';
lfp_EventNames{forcedRedOn} = 'forcedRedOn';
lfp_EventNames{forcedRedOnEM} = 'forcedRedOnEM';
lfp_EventNames{forcedYelOff} = 'forcedYelOff';
lfp_EventNames{forcedYelOffEM} = 'forcedYelOffEM';
lfp_EventNames{forcedYelOn} = 'forcedYelOn';
lfp_EventNames{forcedYelOnEM} = 'forcedYelOnEM';
lfp_EventNames{manualStart} = 'manualStart';
lfp_EventNames{max_reaction_time} = 'max_reaction_time';
lfp_EventNames{postChoCueOn} = 'postChoCueOn';
lfp_EventNames{postFixPtOn} = 'postFixPtOn';
lfp_EventNames{postTargAcq} = 'postTargAcq';
lfp_EventNames{rewardOffAp} = 'rewardOffAp';
lfp_EventNames{rewardOffAv} = 'rewardOffAv';
lfp_EventNames{rewardOnAp} = 'rewardOnAp';
lfp_EventNames{rewardOnAv} = 'rewardOnAv';
lfp_EventNames{squareAcq} = 'squareAcq';
lfp_EventNames{squareOff} = 'squareOff';
lfp_EventNames{targetsOn} = 'targetsOn';
lfp_EventNames{targetsOnEM} = 'targetsOnEM';

lfp_EventColors([airpuffOffAp airpuffOffNoResp AVE_OFF_FORCED]) = {[0.6 0 0]};
lfp_EventColors([airpuffOnAp airpuffOnNoResp AVE_ON_FORCED]) = {[0.8 0 0]};
lfp_EventColors([crossOff squareOff]) = {[0.6 0 0.6]};
lfp_EventColors([rewardOffAp rewardOffAv RWD_OFF_FORCED]) = {[0 0 0.6]};
lfp_EventColors([rewardOnAp rewardOnAv RWD_ON_FORCED]) = {[0 0 0.8]};
lfp_EventColors{CROSS_LEFT} = [0.7 0 0.7];
lfp_EventColors{CROSS_RIGHT} = 'm';
lfp_EventColors{FixPtOn} = 'k';
lfp_EventColors{choiceCueOn} = 'g';
lfp_EventColors{crossAcq} = 'b';
lfp_EventColors{postChoCueOn} = 'c';
lfp_EventColors{postFixPtOn} = 'c';
lfp_EventColors{postTargAcq} = 'c';
lfp_EventColors{squareAcq} = 'r';
lfp_EventColors{targetsOn} = 'k';

lfp_SelectedEventIDs(:) = false;
lfp_SelectedEventIDs([CROSS_LEFT CROSS_RIGHT crossAcq squareAcq ...
    MLTrialStart MLTrialEnd ...
    postChoCueOn postTargAcq postFixPtOn ITI1 ITI2 ...
    choiceCueOn FixPtOn ...
    targetsOn crossOff squareOff ...
    airpuffOnAp airpuffOnNoResp AVE_ON_FORCED ...
    airpuffOffAp airpuffOffNoResp AVE_OFF_FORCED ...
    rewardOnAp rewardOnAv RWD_ON_FORCED ...
    rewardOffAp rewardOffAv RWD_OFF_FORCED ]) = true;


