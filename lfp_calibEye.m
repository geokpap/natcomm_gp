function calibrated = lfp_calibEye(samples, axis, cal)
%calibrated = lfp_calibEye(samples, axis, cal)
% Returns a calibrated eye position waveform, using <samples> as raw
% waveform for <axis> eye position using calibration data <cal>.  <cal> is
% a structure of the type returned by lfp_getEyeCalib.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

oC = 0;       % Cheetah offset
%oC = -47.5;  % Y input measured 8/31/04 with signal generator &
%oscilloscope (but
% then another time it came out +43, so screw that...

%aD = 4096 / 10; % Delphi (behavioral control computer) A/D board gain (32768/2 for new board)
aD = 16384;   % = 32768/2; a priori calc confirmed with signal generator and online 
%                   numeric  display of counts in  Delphi 8/31/04
%aD = 0; % to reveal the location of -g*(oD - offset)
%aD = 5000;   % approx value from peak-to-peak
%aD = 7900;
%aC = 4096 / (4 * 5); % The 5 reflects the recent realization that the output range 
% of eyelink is now one fifth of what it used to be when it was set for the
% old Delphi A/D board
% But wait, that applies equally to the signal going into Delphi A/D, so as
% long as aD, oD are actually factually correct, we shouldn't need to worry
% about the gain of Eyelink's output.
%aC = 4096 / 4; % Cheetah Analog-In (BNC) gain (can also be read from ADBITVOLT from the CSC file-header)
%                 without any attenuators
aC = 790; % Y input measured 8/31/04 with signal generator & oscilloscope @ 4V p-p
%oD = 2048; % Delphi offset (32768 for new board)
oD = 32768;    % considered a priori
%oD = 32400; % crudely smaveraged measured value
%oD = 33500;  % empirically tweaked


% Set offset and gain; gain1 is for negative calibrated coordinate values
% and gain2 for positive values re: offset.  If gain2 is 0, then gain1 is
% applied to all values. Conventions: in raw waveforms, x<0 = right, y<0 =
% up; in calibrated waveforms, x<0 = left, y<0 = down.  Both gains are
% meant to be positive, but no technical issue prevents their being
% negative.

switch axis
    case {'x' 'X'}
        gain1 = cal.XGainL;
        gain2 = cal.XGainR;
        offset = cal.XOffset;
    case {'y' 'Y'}
        gain1 = cal.YGainD;
        gain2 = cal.YGainU;
        offset = cal.YOffset;
end

% Reconstruct the output of the A/D board in the behavioral control
% computer, based on the recorded Cheetah data
if gain2 == 0
    % Use single-valued gain
    calibrated = - gain1 * ...
        ( (samples - oC) * aD/aC + oD - offset );
else
    % Remove offset
    samples = (samples - oC) * aD/aC + oD - offset;
    % Apply proper gain to the proper side of zero
    negvalues = find(samples < 0);
    othervalues = find(samples >= 0);
    calibrated = zeros(size(samples));
    calibrated(negvalues) = - gain1 * samples(negvalues);
    calibrated(othervalues) = - gain2 * samples(othervalues);
end

