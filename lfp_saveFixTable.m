function lfp_saveFixTable(fixations, fname)
% Saves the <fixations2> value returned by lfp_TargSeq to a file named
% <fname> in the same format as lfp_EyeTabulation fixations output file,
% but with an additional column labelled 'targID'.  If <fname> is a
% relative pathname, then the file is saved to the current working
% directory.  Does NOT check for pre-existing <fname>.

%   lfp_SelectedTrials and lfp_BadTrials must be unchanged since
%   <fixations> was generated.  Unfortunately, all that we can do to
%   enforce that condition is check that there is still the same number of
%   trials selected.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

fixoutfid = fopen(fname, 'w');
if fixoutfid == -1
    error('Could not open output file');
end

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= length(fixations)
    error('lfp_saveFixTable:badSelect', ...
        'Trial selection state has changed since generating <fixations>' );
end

fprintf(fixoutfid, 'Trial\tTrialID\tX\tY\tStart\tDuration\tExtra\ttargID\n');
for trialidx = 1:length(trials)
    for row = 1:size(fixations{trialidx},1)
        fprintf(fixoutfid, '%d\t%s\t%d\t%d\t%.6f\t%.6f\t%d\t%d\n', ...
            trials(trialidx), lfp_getTrialID(trials(trialidx)), ...
            fixations{trialidx}(row,:) );
    end
end
fprintf(fixoutfid, '\n');

fclose(fixoutfid);
