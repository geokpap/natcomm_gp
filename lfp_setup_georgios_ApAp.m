% lfp_setup_georgios_ApAp
% Modified from lfp_setup_georgios_ApAv.
%

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
lfp_EventNames{FixPtOnApAp} = 'FixPtOnApAp';
lfp_EventNames{ITI1} = 'ITI1';
lfp_EventNames{ITI2} = 'ITI2';
lfp_EventNames{MLTrialEnd} = 'MLTrialEnd';
lfp_EventNames{MLTrialStart} = 'MLTrialStart';
lfp_EventNames{RWD_OFF_FORCED} = 'RWD_OFF_FORCED';
lfp_EventNames{RWD_ON_FORCED} = 'RWD_ON_FORCED';
lfp_EventNames{TargOffNoResp} = 'TargOffNoResp';
lfp_EventNames{airpuffOnAp} = 'airpuffOnAp';
lfp_EventNames{badTargAcq} = 'badTargAcq';
lfp_EventNames{choiceCueOn} = 'choiceCueOn';
lfp_EventNames{choiceCueOnApAp} = 'choiceCueOnApAp';
lfp_EventNames{crossAcq} = 'crossAcq';
lfp_EventNames{crossOff} = 'crossOff';
lfp_EventNames{crossAcqApAp} = 'crossAcqApAp';
lfp_EventNames{event_zero} = 'event_zero';
lfp_EventNames{fixbreak} = 'fixbreak';
lfp_EventNames{fixptOffFixbreak} = 'fixptOffFixbreak';
lfp_EventNames{forcedRedOff} = 'forcedRedOff';
lfp_EventNames{forcedRedOn} = 'forcedRedOn';
lfp_EventNames{forcedYelOff} = 'forcedYelOff';
lfp_EventNames{forcedYelOn} = 'forcedYelOn';
lfp_EventNames{manualStart} = 'manualStart';
lfp_EventNames{max_reaction_time} = 'max_reaction_time';
lfp_EventNames{postChoCueOn} = 'postChoCueOn';
lfp_EventNames{postFixPtOn} = 'postFixPtOn';
lfp_EventNames{postTargAcq} = 'postTargAcq';
lfp_EventNames{rewardOffAv} = 'rewardOffAv';
lfp_EventNames{rewardOffRedApAp} = 'rewardOffRedApAp';
lfp_EventNames{rewardOffYelApAp} = 'rewardOffYelApAp';
lfp_EventNames{rewardOnAv} = 'rewardOnAv';
lfp_EventNames{rewardOnRedApAp} = 'rewardOnRedApAp';
lfp_EventNames{rewardOnYelApAp} = 'rewardOnYelApAp';
lfp_EventNames{squareOff} = 'squareOff';
lfp_EventNames{squareAcqYelApAp} = 'squareAcqYelApAp';
lfp_EventNames{targetsOffYelApAp} = 'targetsOffYelApAp';
lfp_EventNames{targetsOn} = 'targetsOn';

lfp_EventColors([crossOff squareOff]) = {[0.6 0 0.6]};
lfp_EventColors{CROSS_LEFT} = [0.7 0 0.7];
lfp_EventColors{CROSS_RIGHT} = 'm';
lfp_EventColors{FixPtOn} = 'k';
lfp_EventColors{FixPtOnApAp} = 'k';
lfp_EventColors{choiceCueOn} = 'g';
lfp_EventColors{choiceCueOnApAp} = 'g';
lfp_EventColors{crossAcq} = 'b';
lfp_EventColors{crossAcqApAp} = 'b';
lfp_EventColors{postChoCueOn} = 'c';
lfp_EventColors{postFixPtOn} = 'c';
lfp_EventColors{postTargAcq} = 'c';
lfp_EventColors{squareAcqYelApAp} = 'r';
lfp_EventColors{targetsOn} = 'k';



