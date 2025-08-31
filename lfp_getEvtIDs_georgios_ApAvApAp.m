% lfp_getEvtIDs_georgios_ApAvApAp
%NOTE
% Due to a conflict in the MonkeyLogic code over whether event 213 means
% FixPtOnApAp or RWD_ON_FORCED, I replace all forced rwd-on events with evt
% 210 in 'lfp_otherPreprocess' regardless of what block they are in and
% what cue they refer to, and abolish what I was going to call
% RWD_ON_FORCEDYelApAp (originally 210) and RWD_ON_FORCEDRedApAp
% (originally 250).

%$Rev:  $
%$Date:  $
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
CROSS_LEFT = 251;
CROSS_RIGHT = 252;
FORCED = 254;
FixPtOn = 214;
FixPtOnApAp = 213;
MLTrialEnd = 18;
MLTrialStart = 9;
RWD_OFF_FORCED = 227;
RWD_OFF_FORCEDYelApAp = 217;
RWD_OFF_FORCEDRedApAp = 250;
RWD_ON_FORCED = 210;
TargOffNoResp = 203;

airpuffOffAp = 240;
airpuffOffNoResp = 236;
airpuffOnAp = 239;
airpuffOnNoResp = 235;
badTargAcq = 206;
choiceCueOn = 233;
choiceCueOnApAp = 224;
crossAcq = 219;
crossAcqApAp = 218;
crossOff = 201;
crossOffRedApAp = 226;
event_zero = 666;
fixbreak = 204;
fixptOffFixbreak = 205;
forcedRedOff = 216;
forcedRedOffApAp = 209;
forcedRedOn = 215;
forcedRedOnApAp = 208;
forcedYelOff = 212;
forcedYelOffApAp = 223;
forcedYelOn = 211;
forcedYelOnApAp = 222;
manualStart = 255;
max_reaction_time = 202;
rewardOffRedApAp = 221;
rewardOffYelApAp = 247;
rewardOffAv = 245;
rewardOffAp = 242;
rewardOnRedApAp = 220;
rewardOnYelApAp = 253;
rewardOnAv = 244;
rewardOnAp = 241;
squareAcq = 243;
squareAcqYelApAp = 237;
squareOff = 207;
targetsOn = 234;
targetsOnApAp = 225;
targetsOffYelApAp = 248;

% Events to calculate:
postChoCueOn = 256;
postTargAcq = 257;
postFixPtOn = 258;
ITI1 = 259;
ITI2 = 260;

