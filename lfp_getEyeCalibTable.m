function table = lfp_getEyeCalibTable(option)
%table = lfp_getEyeCalibTable
%   If lfp_ManualEyeCalib is not empty, return its value (note that
%   this value might not be sorted).  Otherwise, construct a table in the
%   same format as lfp_ManualEyeCalib containing the calibration values
%   extracted from lfp_TrialParams for all trials.  Duplicate calibrations
%   are omitted, i.e. the table returned contains only the changes in the
%   calibrations.  In all cases the table is returned sorted by timestamp.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

forceflag = false;
if nargin > 0
    switch option
        case 'force'
            forceflag = true;
        otherwise
            error('lfp_getEyeCalibTable:unknownOpt', ...
                'The option %s is not recognized', option );
    end
end

if ~forceflag && ~isempty(lfp_ManualEyeCalib)
    table = lfp_ManualEyeCalib;
    disp('Using eye calibrations in lfp_ManualEyeCalib');
else
    table = {};
    for trial = 1:size(lfp_TrialIndex,1)
        if forceflag
            c = lfp_getEyeCalib(trial, 'force');
        else
            c = lfp_getEyeCalib(trial);
        end
        if isempty(table) || ~isequal(c, table{end,2})
            table{end+1,2} = c;
            table{end,1} = lfp_Events(lfp_TrialIndex(trial,1), 1);
        end
    end
    disp('Using eye calibrations from trial parameters');
end
keycol = cell2mat(table(:,1));
[keys, index] = sortrows(keycol);
table = table(index,:);
