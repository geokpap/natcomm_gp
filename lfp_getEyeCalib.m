function cal = lfp_getEyeCalib(trial, varargin)
%cal = lfp_getEyeCalib(trial)
%cal = lfp_getEyeCalib(..., 'force')
%cal = lfp_getEyeCalib(..., 'v')
% Returns a structure with elements
%   XOffset
%   XGainL
%   XGainR
%   YOffset
%   YGainU
%   YGainD
% XGainL and YGainD are the fields used for symmetric gain values.  If
% lfp_ManualEyeCalib contains values, the last calibration before the
% current trial is returned; if there is none, then the first calibration
% is returned.  If lfp_ManualEyeCalib is empty, the values are retrieved
% from the trial parameters (aka BD section).
% OPTIONS :
%   'force'
% The calibrations are returned from the trial parameters even if there are
% values in lfp_ManualEyeCalib.
%   'v'
% Verbose, displays trial number and calibration.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

forceflag = false;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'char')
        switch varargin{argnum}
            case 'force'
                forceflag = true;
            case 'v'
                verboseflag = true;
            otherwise
                error('lfp_getEyeCalib:unknownOpt', ...
                    'The option "%s" is not recognized', option );
        end
    end
    argnum = argnum + 1;
end

if isempty(lfp_ManualEyeCalib) || forceflag
    % extract vals from BDoutput for the current trial
    if lfp_TrialParams{trial}(1) < 7
        error('lfp_getEyeCalib:unknownFormat1', ...
            'Unknown format: %d', lfp_TrialParams{trial}(1) );
    elseif lfp_TrialParams{trial}(1) == 7   % old Joey format
        cal.XOffset = lfp_TrialParams{trial}(21) + ...
            lfp_TrialParams{trial}(22) * 256;
        cal.XGainL = lfp_TrialParams{trial}(29) / ...
            lfp_TrialParams{trial}(30);
        cal.XGainR = 0;
        
        cal.YOffset = lfp_TrialParams{trial}(23) + ...
            lfp_TrialParams{trial}(24) * 256;
        cal.YGainD = lfp_TrialParams{trial}(31) / ...
            lfp_TrialParams{trial}(32);
        cal.YGainU = 0;
    else
        % 2-byte format code.
        % Construct format code from this and the next param, treating
        % them respectively as the high-order and low-order pairs of hex
        % digits of a single number.
        formatnum = 256*bitand(255, lfp_TrialParams{trial}(1)) + ...
            bitand(255, lfp_TrialParams{trial}(2));
        switch formatnum
            case 2817 % 0xB01: start new Joey formats; 
                % corrected gain formulae to match old Delphi code 9/11/2004
                cal.XOffset = lfp_TrialParams{trial}(22) + ...
                    lfp_TrialParams{trial}(23) * 256;
                cal.XGainL = 0.5 * lfp_TrialParams{trial}(30) / ...
                    lfp_TrialParams{trial}(31);
                cal.XGainR = 0;
                
                cal.YOffset = lfp_TrialParams{trial}(24) + ...
                    lfp_TrialParams{trial}(25) * 256;
                cal.YGainD = 0.5 * lfp_TrialParams{trial}(32) / ...
                    lfp_TrialParams{trial}(33);
                cal.YGainU = 0;
                
            case {2818, 2819, 2820, 2821, 2822, 2823, 2824, 2825, 2826,...
                   2827, 2828, 2829, 2830, 2831}
                % 0xB02 & on (new formats for corrected Delphi
                % display code starting 9/10/2004)
                cal.XOffset = lfp_TrialParams{trial}(22) + ...
                    lfp_TrialParams{trial}(23) * 256;
                cal.XGainL = lfp_TrialParams{trial}(30) / ...
                    lfp_TrialParams{trial}(31);
                cal.XGainR = 0;

                cal.YOffset = lfp_TrialParams{trial}(24) + ...
                    lfp_TrialParams{trial}(25) * 256;
                cal.YGainD = lfp_TrialParams{trial}(32) / ...
                    lfp_TrialParams{trial}(33);
                cal.YGainU = 0;
                
            case 2309   % 0x905
                XCalibEye = (lfp_TrialParams{trial}(38)-hex2dec('2600'))*256 + ...
                    (lfp_TrialParams{trial}(39)-hex2dec('2700'));
                XGainEye1L = lfp_TrialParams{trial}(42)-hex2dec('2A00');
                XGainEye1R = lfp_TrialParams{trial}(44)-hex2dec('2C00');
                XGainEye2L = lfp_TrialParams{trial}(43)-hex2dec('2B00');
                XGainEye2R = lfp_TrialParams{trial}(45)-hex2dec('2D00');
                YCalibEye = (lfp_TrialParams{trial}(40)-hex2dec('2800'))*256 ...
                    + (lfp_TrialParams{trial}(41)-hex2dec('2900'));
                YGainEye1D = lfp_TrialParams{trial}(48)-hex2dec('3000');
                YGainEye1U = lfp_TrialParams{trial}(46)-hex2dec('2E00');
                YGainEye2D = lfp_TrialParams{trial}(49)-hex2dec('3100');
                YGainEye2U = lfp_TrialParams{trial}(47)-hex2dec('2F00');
                cal.XOffset = 65535 - XCalibEye;
                cal.YOffset = YCalibEye;
                cal.XGainL = XGainEye1L / XGainEye2L;
                cal.XGainR = XGainEye1R / XGainEye2R;
                cal.YGainU = YGainEye1U / YGainEye2U;
                cal.YGainD = YGainEye1D / YGainEye2D;

            case 2310   % 0x906
                XCalibEye = (lfp_TrialParams{trial}(38)-hex2dec('2600'))*256 + ...
                    (lfp_TrialParams{trial}(39)-hex2dec('2700'));
                XGainEye1L = lfp_TrialParams{trial}(42)-hex2dec('2A00');
                XGainEye1R = lfp_TrialParams{trial}(44)-hex2dec('2C00');
                XGainEye2L = lfp_TrialParams{trial}(43)-hex2dec('2B00');
                XGainEye2R = lfp_TrialParams{trial}(45)-hex2dec('2D00');
                YCalibEye = (lfp_TrialParams{trial}(40)-hex2dec('2800'))*256 ...
                    + (lfp_TrialParams{trial}(41)-hex2dec('2900'));
                YGainEye1D = lfp_TrialParams{trial}(48)-hex2dec('3000');
                YGainEye1U = lfp_TrialParams{trial}(46)-hex2dec('2E00');
                YGainEye2D = lfp_TrialParams{trial}(49)-hex2dec('3100');
                YGainEye2U = lfp_TrialParams{trial}(47)-hex2dec('2F00');
                cal.XOffset = 65535 - XCalibEye;
                cal.YOffset = YCalibEye;
                cal.XGainL = XGainEye1L / XGainEye2L;
                cal.XGainR = XGainEye1R / XGainEye2R;
                cal.YGainU = YGainEye1U / YGainEye2U;
                cal.YGainD = YGainEye1D / YGainEye2D;
                
            case 2311   % 0x907
                XCalibEye = (lfp_TrialParams{trial}(38)-hex2dec('2600'))*256 + ...
                    (lfp_TrialParams{trial}(39)-hex2dec('2700'));
                XGainEye1L = lfp_TrialParams{trial}(42)-hex2dec('2A00');
                XGainEye1R = lfp_TrialParams{trial}(44)-hex2dec('2C00');
                XGainEye2L = lfp_TrialParams{trial}(43)-hex2dec('2B00');
                XGainEye2R = lfp_TrialParams{trial}(45)-hex2dec('2D00');
                YCalibEye = (lfp_TrialParams{trial}(40)-hex2dec('2800'))*256 ...
                    + (lfp_TrialParams{trial}(41)-hex2dec('2900'));
                YGainEye1D = lfp_TrialParams{trial}(48)-hex2dec('3000');
                YGainEye1U = lfp_TrialParams{trial}(46)-hex2dec('2E00');
                YGainEye2D = lfp_TrialParams{trial}(49)-hex2dec('3100');
                YGainEye2U = lfp_TrialParams{trial}(47)-hex2dec('2F00');
                cal.XOffset = 65535 - XCalibEye;
                cal.YOffset = YCalibEye;
                cal.XGainL = XGainEye1L / XGainEye2L;
                cal.XGainR = XGainEye1R / XGainEye2R;
                cal.YGainU = YGainEye1U / YGainEye2U;
                cal.YGainD = YGainEye1D / YGainEye2D;

            case {2312, 2313, 2314, 2315, 2316, 2317, 2318}   % 0x908 to 0x90E
                % At the time of this rewrite, this was the last format
                % number, so I have rewritten this one in shorter form.
                % (Or so I thought...)
                %   DG 12-Jan-2005
                % TMD updated format numbers 9/21/07
                cal.XOffset = 65535 - (bitand(lfp_TrialParams{trial}(38), 255)*256 + ...
                    bitand(lfp_TrialParams{trial}(39), 255));
                cal.YOffset = bitand(lfp_TrialParams{trial}(40), 255)*256 ...
                    + bitand(lfp_TrialParams{trial}(41), 255);
                cal.XGainL = bitand(lfp_TrialParams{trial}(42), 255) / ...
                    bitand(lfp_TrialParams{trial}(43), 255);
                cal.XGainR = bitand(lfp_TrialParams{trial}(44), 255) / ...
                    bitand(lfp_TrialParams{trial}(45), 255);
                cal.YGainU = bitand(lfp_TrialParams{trial}(46), 255) / ...
                    bitand(lfp_TrialParams{trial}(47), 255);
                cal.YGainD = bitand(lfp_TrialParams{trial}(48), 255) / ...
                    bitand(lfp_TrialParams{trial}(49), 255);
                
                    
            otherwise
                error('lfp_getEyeCalib:unknownFormat2', ...
                    'Unknown format: 0x%X', formatnum );
        end
    end
    
else
    % Use values saved in lfp_ManualEyeCalib (which should be
    % automatically loaded by lfp_read from lfp_ManualEyeCalib.Mat).
    % Find the last calibration before the current trial:
    trialstarttime = lfp_Events(lfp_TrialIndex(trial,1),1);
    eyecalibtimes = cell2mat(lfp_ManualEyeCalib(:,1));
    TimeDifference = trialstarttime - eyecalibtimes;
    if all(TimeDifference < 0)
        % There is no calibration before the current trial; use the
        % first calibration.
        calibrationTSindex = 1;
    else
        minTD = min(TimeDifference(find(TimeDifference >= 0)));
        calibrationTSindex = find(TimeDifference == minTD);
    end
    
    cal.XOffset = lfp_ManualEyeCalib{calibrationTSindex,2}.XOffset;
    cal.XGainL = lfp_ManualEyeCalib{calibrationTSindex,2}.XGainL;
    cal.XGainR = lfp_ManualEyeCalib{calibrationTSindex,2}.XGainR;
    
    cal.YOffset = lfp_ManualEyeCalib{calibrationTSindex,2}.YOffset;
    cal.YGainU = lfp_ManualEyeCalib{calibrationTSindex,2}.YGainU;
    cal.YGainD = lfp_ManualEyeCalib{calibrationTSindex,2}.YGainD;
end
if verboseflag
    if forceflag
        disp('Forced reading eye calib from trial params')
    end
    disp(sprintf('Eye calibs for trial %d:', trial));
    disp(cal);
end

