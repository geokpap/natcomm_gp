function [seqs, fixations2] = lfp_TargSeq(fixations, targets, ...
    startEvts, stopEvts, RT, varargin)
%LFP_TARGSEQ extracts sequence of targets fixated for each trial
%seqs = lfp_TargSeq(fixations, targets, startEvts, stopEvts, RT)
%seqs = lfp_TargSeq(..., 'inprogress')
%seqs = lfp_TargSeq(..., 'startoffset', offset)
%seqs = lfp_TargSeq(..., 'stopoffset', offset)
%   lfp_SelectedTrials and lfp_BadTrials must be unchanged since
%   <fixations> was generated.  Unfortunately, all that we can do to
%   enforce that condition is check that there is still the same number of
%   trials selected.

%seqs = lfp_TargSeq(fixations, targets, startEvts, stopEvts, RT)
%   For each element of <fixations>, extracts sequence of targets fixated.
%   <fixations> is in the format returned by lfp_EyeTabulation. The eye
%   coordinates are assumed to be calibrated. <targets> contains target
%   center X and Y in the first two columns and target radius in the third
%   for an arbitrary number of targets (one row per target).  It is assumed
%   that the target radii are small enough so that targets do not overlap.
%   <startEvts> and <stopEvts> define the time interval over which the
%   target sequence is extracted; they are lists of event IDs, and the
%   first trial event that is a member of the list is used.  We ignore any
%   fixations that begin earlier than the time of start event plus <RT>
%   (Reaction Time), and include the fixation-in-progress at time of stop
%   event.  <seqs> is a cell vector containing a numeric array for each
%   trial.  Each numeric array contains the fixation start times in the
%   first column; the second column contains integers representing the
%   row of <targets> for the target fixated; the third column contains
%   fixation end times.  If the fixation is not in any target, then the
%   value is 0.  Therefore, a cell in <seqs> may contain [] only if there
%   are no fixations in the trial between the start and stop events
%   (unlikely).
%[seqs, fixations2] = ...
%       lfp_TargSeq(fixations, targets, startEvts, stopEvts, RT)
%   <fixations2> is a modified version of <fixations> that has target IDs
%   in an additional column at the end, with NaN for any target that is not
%   included in <seqs>.  In the case of a trial where no fixations are
%   included, fixations2{trial} is NOT empty, but contains a single NaN.
%OPTIONS
%   'inprogress' - include the fixation-in-progress at the time of the
%       start event (this makes RT irrelevant, but it must still be given a
%       value to keep the correct argument count)
%   'startoffset' - must be followed by an offset in msec that is added to
%       the timestamp of the start event to set the beginning of the time
%       interval (see <startEvts>).  The offset may be negative.
%   'stopoffset' - like 'startoffset', but applies to stop event.

% 10/19/04 DG changed errors to warning...continue in 65-75.
% 7/25/05 DG added 3rd column to output

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 5
    error('lfp_TargSeq:missingArg', ...
        'Insufficient arguments; type "help lfp_TargSeq".' );
end

fixations2 = fixations;

argnum = 1;
inprogflag = false;
startoffset = 0;
stopoffset = 0;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'inprogress'
            inprogflag = true;
        case 'startoffset'
            startoffset = 1e-3*varargin{argnum+1};
            argnum = argnum+1;
        case 'stopoffset'
            stopoffset = 1e-3*varargin{argnum+1};
            argnum = argnum+1;
        otherwise
            error('lfp_fragmentFiles:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end


trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= length(fixations)
    error('lfp_TargSeq:badSelect', ...
        'Trial selection state has changed since generating <fixations>' );
end

for trialidx = 1:length(fixations)
    seqs{trialidx} = [];
    fixations2{trialidx}(:, end+1) = NaN;
    evtrange = lfp_TrialIndex(trials(trialidx), 1) : lfp_TrialIndex(trials(trialidx), 2);
    startTS = startoffset + lfp_Events(find(...
        ismember(lfp_Events(evtrange, 2), startEvts) ) + evtrange(1) - 1, 1);
    stopTS = stopoffset + lfp_Events(find(...
        ismember(lfp_Events(evtrange, 2), stopEvts) ) + evtrange(1) - 1, 1);
    if isempty(startTS)
        warning('lfp_TargSeq:nostart', ...
            'No start event in trial %d', trials(trialidx) );
        continue
    end
    if isempty(stopTS)
        warning('lfp_TargSeq:nostop', ...
            'No stop event in trial %d', trials(trialidx) );
        continue
    end
    for fixidx = 1:size(fixations{trialidx}, 1)
        fixTS = fixations{trialidx}(fixidx,3);
        fixend = fixTS + fixations{trialidx}(fixidx,4);
        if inprogflag
            include = (fixend > startTS) && (fixTS < stopTS);
        else
            include = (fixTS > startTS + RT) && (fixTS < stopTS);
        end
        if include
            x = fixations{trialidx}(fixidx,1);
            y = fixations{trialidx}(fixidx,2);
            targID = lfp_xy2targID(x, y, targets);
            seqs{trialidx}(end+1,:) = [fixTS targID fixend];
            fixations2{trialidx}(fixidx,end) = targID;
        elseif fixTS >= stopTS
            % remaining fixations are also after stopTS
            break
        end
    end
end
