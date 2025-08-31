function saveSeqTables(sessiondir, varargin)
% saveSeqTables(sessiondir)
% saveSeqTables(..., 'fixations', fixations)

% Accept fixation table as optional argument,
% or look for it in the sessiondir, or as a last resort recalculate it.

% Color sequences not implemented yet.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

fixationsflag = false;
saccadesflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'fixations'
            fixationsflag = true;
            argnum = argnum + 1;
            fixations = varargin{argnum};
        otherwise
            error('saveSeqTables:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[pathstr,fragname,ext,versn] = fileparts(sessiondir);
lfp_read2('preset', sessiondir, ...
    {'events.EVTSAV', 'calsdgX.MAT', 'calsdgy.MAT', 'v.MAT'} );
lfp_getEvtIDs;

if ~fixationsflag
    % load or calculate fixations table
    fixationsfile = fullfile(pathstr, fragname, 'fixsacc.mat');
    if exist(fixationsfile) == 2
        load(fixationsfile);
    else
        [fixations, saccades] = lfp_EyeTabulation(1,2,3);
        save(fullfile(sessiondir, 'fixsacc.mat'), 'fixations', 'saccades');
    end
end
    
% should change the following section of code to create multiple pairs of
% raw and rank-ordered sequences between relevant behavioral events, not
% merely between H1Start and M1out.  In order to do that, we must change
% both the immediate line as well as the call to lfp_SeqTabulation, since
% it presently only saves rawseq and s4seq.

% Each row of <boundses> contains a start event and end event:
boundses = [ H0Start H1Start
    H1Start M1out
    M1out H2Start
    H2Start M2out
    M2out H3Start
    H3Start M3out
    M3out H4Start
    H4Start PRDelay
    PRDelay RWD
    RWD ITI ];
intervalnames = { 'H0'
    'H1'
    'M1'
    'H2'
    'M2'
    'H3'
    'M3'
    'H4'
    'PRD'
    'RWD' };

for intervalnum = 1:length(intervalnames)
    % NOTE: the following call to lfp_TargSeq determines the target of each
    % fixation based on the target information supplied in lfp_targets_joey,
    % which is only valid from 6/2004
    seqs(:,intervalnum) = reshape(...
        lfp_TargSeq(fixations, lfp_targets_joey, ...
        boundses(intervalnum,1), boundses(intervalnum,2), ...
        0, 'inprogress'), ...
        [], 1);
    colnames{(intervalnum-1)*4 + 1} = ...
        [intervalnames{intervalnum} 'i'];
    colnames{(intervalnum-1)*4 + 2} = ...
        [intervalnames{intervalnum} 'r0'];
    colnames{(intervalnum-1)*4 + 3} = ...
        [intervalnames{intervalnum} 'r'];
    colnames{(intervalnum-1)*4 + 4} = ...
        [intervalnames{intervalnum} 'PT'];
end

% save the data
lfp_SeqTabulation(seqs, 'headers', colnames, 'renumno0', 'preset', ...
    fullfile(sessiondir, 'seq.txt') );
lfp_TrialParamsTabulation('preset', ...
    fullfile(sessiondir, 'tp.txt') );