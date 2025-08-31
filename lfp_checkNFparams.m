function badtrials = lfp_checkNFparams(params)

%$Rev: 111 $
%$Date: 2010-02-22 00:21:09 -0500 (Mon, 22 Feb 2010) $
%$Author: dgibson $

lfp_declareGlobals;

tp = cell2mat(params');
isbad = false(size(tp));
isbad(:, 1:6) = tp(:, 1:6) == -1;
isbad(:, [7 8 12]) = ~ismember(tp(:, [7 8 12]), [0 1]);
isbad(:, 9) =  tp(:, 9) == -1;
isbad(:, 10) = ~ismember(tp(:, 10), [0:7 14:17 19:26]); % from getConditionDescription.m
isbad(:, 11) = ~ismember(tp(:, 11), 0:7);
badtrials = find(any(isbad,2));
if ~isempty(badtrials)
    report = sprintf('Session %s contains bad trial params\nbad trials / total = %d / %d = %d%%\n', ...
        lfp_SessionNames{1}, length(badtrials), ...
        length(lfp_SelectedTrials), ...
        round(100*length(badtrials)/length(lfp_SelectedTrials)));
    for trial = reshape(badtrials, 1, [])
        report = sprintf('%strial %d bad params %s\n', ...
            report, trial, mat2str(find(isbad(trial, :))) );
    end
    lfp_log(report);
    if length(badtrials) > 50
        warnmsg = 'More than 50 bad trials; see lfp_lib.log';
    else
        warnmsg = report;
    end
    warning('lfp_checkNFparams:bad', '%s', warnmsg);
end
