function [fixO, fixT, fixY, fixR, fixG, fixB] = lfp_fixtrains
% generates a histogram of fixation "counts" for each of the possible 6 of
% fixation targets

% NOTE: change the function to accept a fix table if it is provided, and to
% compute it from scratch if it isn't

% Written by Joey Feingold spring 2005
% DG modified 7/25/05 for new lfp_TargSeq output format

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

lfp_getEvtIDs;

warning off lfp_EyeTabulation:nostop
%generate fix and sacc tables for all selected trials (will place the
%entire, unedited, fix and sacc list for any trial lacking either the start
%or the stop event):
[fix, sacc] = lfp_EyeTabulation(1,2,3);

warning off lfp_TargSeq:nostop
% trial selection must remain constant until the following line, but may
% change thereafter. lfp_TargSeq will return empty for trials missing either
% the start or stop event
seqs = lfp_TargSeq(fix, lfp_targets_joey, H0Start, BDOutputStart, 0, 'inprogress');

% Now generate sequences of the ordered (R-11, G-12, B-13, Y-9) targets
cueseqs = lfp_RenumberTargSeq(seqs, 10);

resolution = .01;

fixO = []; % all "other" fixations
fixY = []; % fixation of Yellow cue
fixR = []; % fixation of Red cue
fixG = []; % fixation of Green cue
fixB = []; % fixation of Blue cue
fixT = []; % fixation of a non-target (irrelevant) cue

fixtargetIDs = [0 1 2 3 4 5 6 7 8 9 11 12 13];
finalfix = cell(length(fixtargetIDs),1);

for i = 1:length(cueseqs)
    for j = 1:length(fixtargetIDs)
        fixOnOff = cueseqs{i}(find(cueseqs{i}(:,2) == fixtargetIDs(j)),[1 3]);
        for k = 1:size(fixOnOff,1)
            fixexpanded = [fixOnOff(k,1) : resolution : fixOnOff(k,2)];
            finalfix{j} = [finalfix{j} fixexpanded];
        end
    end
end

fixO = finalfix{1};

for i = 2:9
    fixT = [fixT finalfix{i}]; 
end    

fixY = finalfix{10};
fixR = finalfix{11};
fixG = finalfix{12};
fixB = finalfix{13};

save(fullfile(lfp_DataDir, 'fixO'), 'fixO');
save(fullfile(lfp_DataDir, 'fixT'), 'fixT');
save(fullfile(lfp_DataDir, 'fixY'), 'fixY');
save(fullfile(lfp_DataDir, 'fixR'), 'fixR');
save(fullfile(lfp_DataDir, 'fixG'), 'fixG');
save(fullfile(lfp_DataDir, 'fixB'), 'fixB');