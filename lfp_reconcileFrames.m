function canonicalstamps = lfp_reconcileFrames(timestamps)
% When Cheetah recording stops, it sometimes happens that some channels
% record one frame more than the other channels.  Rodent data may contain
% many starts and stops, so there can be frame discrepancies in the middles
% of files.  We make the files uniform by deleting the extra frames
% wherever they occur.  Returns the list of timestamps common to all files.
% Note: this function requires that lfp_Samples be in "sample x frame"
% format.

%$Rev: 81 $
%$Date: 2009-09-03 17:52:27 -0400 (Thu, 03 Sep 2009) $
%$Author: dgibson $

lfp_declareGlobals;

[canonicalstamps, activefnums] = dg_reconcileFrames(timestamps);

if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', 'Reconciling Frames');
    WaitBarSteps = 2*length(activefnums);
    CurrentWaitBarStep = 0;

    CurrentWaitBarStep = CurrentWaitBarStep + 1;
    waitbar(CurrentWaitBarStep/WaitBarSteps, hWaitBar);
end
for filenum = activefnums
    if ~isempty(timestamps{filenum})
        mismatches = true;
        while ~isempty(mismatches)
            lengthdiff = length(timestamps{filenum}) - ...
                length(canonicalstamps); % always >= 0
            mismatches = find(timestamps{filenum}(1:end-lengthdiff) ...
                ~= canonicalstamps );
            if ~isempty(mismatches)
                timestamps{filenum}(mismatches(1)) = [];
                lfp_Samples{filenum}(:,mismatches(1)) = [];
            elseif lengthdiff
                timestamps{filenum}(end-lengthdiff+1:end) = [];
                lfp_Samples{filenum}(:,end-lengthdiff+1:end) = [];
            end
        end
        if ~lfp_NoWaitbar
            CurrentWaitBarStep = CurrentWaitBarStep + 1;
            waitbar(CurrentWaitBarStep/WaitBarSteps, hWaitBar);
        end
    end
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end
