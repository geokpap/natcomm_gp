% lfp_getJoyCalibs
% Assumes that lfp_lib global variables are accessible and that the
% variable <trial> has a valid value for a trial number.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

switch lfp_TrialParams{trial}(1)
    case 7
        XCalib = lfp_TrialParams{trial}(25) + lfp_TrialParams{trial}(26) * 256;
        YCalib = lfp_TrialParams{trial}(27) + lfp_TrialParams{trial}(28) * 256;
    otherwise
        error('lfp_getJoyCalib:unknownFormat', ...
            'Unknown format: %d', lfp_TrialParams{trial}(1) );
end
