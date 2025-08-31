%test_lfp_burstAnalysis
% Thia was hanging around in my MATLAB directory with a modification date
% of Apr 20  2012, so I thought I should check it in.  -DG

%$Rev: 305 $
%$Date: 2013-09-05 18:24:19 -0400 (Thu, 05 Sep 2013) $
%$Author: dgibson $

if matlabpool('size') == 0
    matlabpool;
end
lfp_declareGlobals
lfp_changeSetup('joey_H0_1000');
lfp_getEvtIDs;

fragdir = '/smbshare/tmp/dgibson/testdata/hh092807/controls'; 
[implant_str, group_tt_alt_flag] = getImplantstr(upper('hh092807'));
cscfn = 'e30.mat';
evtsfname = 'eventsh0-1000_epochs.evtsav';
task = 'ZeroM3T'; 
rule = 'all_correct'; 
lfp_read2('preset', fragdir, {evtsfname});
lfp_Events(lfp_Events(:,2)==124, 2) = 127;
myrule = setRule_New(rule, task, group_tt_alt_flag); 
lfp_selectByRule(myrule, sprintf('%s %s', rule, task)); 
sum(lfp_SelectedTrials)

disp('REMEMBER TO SET aligns!')
% Sing:
aligns = {H0Start H1Start 991 992 RWD ITIon};
% Sim:
aligns = {H0Start H1Start 991 H2Start M2out H3Start ...
    M3out 992 RWD ITIon};
% ZeroM3T:
aligns = {H0Start H1Start RWD ITIon};
% OneM3T:
aligns = {H0Start H1Start 991 992 RWD ITIon};

lfp_XLimAll = [-1 3]; lfp_AlignmentRef = aligns{1}; lfp_selectByDuration('and');
sum(lfp_SelectedTrials)

[values, plotdata] = lfp_makepasteup(aligns, ...
@lfp_burstAnalysis, [], lfp_ActiveFilenums(2:end), [], ...
'binsize', 10, 'keepboot', ...
'boot', 1000, 'numsets', 1, 'numshuffles', 1000, ...
'plevel', .025, 'clevel', .025, ...
'pasteup_opts', {'data', 'autoselect'});

disp('SAVE values AND plotdata!')
save(sprintf('/smbshare/tmp/dgibson/testdata/coact/%s_%s.mat', rule, ...
task), 'values', 'plotdata')