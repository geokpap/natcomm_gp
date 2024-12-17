% lfp_getEvtIDs_georgios_ApAv_rev_only
% The brute force "make it work" event codes, including codes from
% "Eventmarkers_ApAv.txt".  Works with session types 'ApAv' and 'ApAv_em'.

%$Rev: 45 $
%$Date: 2023-03-31 20:14:39 -0400 (Fri, 31 Mar 2023) $
%$Author: dgibson $

CSCFileRegexp = '.*down\d+\.';
CSCFileExt = 'mat';
EVFilename1 = 'events.nev';
EVFilename2 = 'just say no:\\\.';
useFileSelect = true;
read_mode = 'mat';
usestrobe=false;

AVE_OFF_FORCED = 228;
AVE_ON_FORCED = 238;
CROSS_LEFT = 251; % checked
CROSS_RIGHT = 252; % checked
FORCED = 254; % checked
FixPtOn = 230; % checked
MLTrialEnd = 18; % checked
MLTrialStart = 9; % checked
RWD_OFF_FORCED = 227;
RWD_ON_FORCED = 217;
TargOffNoResp = 203;
airpuffOffAp = 240;
airpuffOffNoResp = 236;
airpuffOnAp = 239;
airpuffOnNoResp = 235;
badTargAcq = 206;
choiceCueOn = 233;
choiceCueOnEM = 246;
crossAcq = 219;
crossOff = 201;
event_zero = 666;
fixbreak = 204;
fixptOffFixbreak = 205;
forcedRedOff = 216;
forcedRedOffEM = 230;
forcedRedOn = 215;
forcedRedOnEM = 229;
forcedYelOff = 212;
forcedYelOffEM = 232;
forcedYelOn = 211;
forcedYelOnEM = 231;
manualStart = 255;
max_reaction_time = 202;
rewardOffAp = 242;
rewardOffAv = 245;
rewardOffRedApAp = 221;
rewardOffYelApAp = 247;
rewardOnAp = 241;
rewardOnAv = 244;
rewardOnRedApAp = 220;
rewardOnYelApAp = 253;
squareAcq = 243;
squareOff = 207;
targetsOn = 234;
targetsOnEM = 247;

% Events to calculate:
postChoCueOn = 256;
postTargAcq = 257;
postFixPtOn = 258;
ITI1 = 259;
ITI2 = 260;

