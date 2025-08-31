function lfp_markDataWindows(windows, ch)
% Replaces the CSC wave in channel number <ch> with a new one containing
% ones at all data points included in <windows>, zeros elsewhere.  If <ch>
% is empty or not given, then creates a new channel.

%$Rev: 75 $
%$Date: 2009-07-21 13:55:56 -0400 (Tue, 21 Jul 2009) $
%$Author: dgibson $

if nargin < 2
    ch = [];
end

lfp_declareGlobals;

x = zeros(size(lfp_Samples{lfp_ActiveFilenums(1)}));
for row = 1:size(windows,1)
    x(windows(row,1):windows(row,2)) = 1;
end
if isempty(ch)
    lfp_createWave(@lfp_waverecord, lfp_ActiveFilenums(1), x, ...
        'name', 'winsamples');
else
    lfp_createWave(@lfp_waverecord, lfp_ActiveFilenums(1), x, ...
        'name', 'winsamples', ...
        'replace', ch);
end
