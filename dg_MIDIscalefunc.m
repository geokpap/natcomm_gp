function values = dg_MIDIscalefunc(values)

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

% b is a scale param that determines "how logarithmic" the log part is.

zp = 60;    % zero point, i.e. the output produced for 0 in.
gain = 6;
xjoint = 1; % the input value at which we switch from linear to log

% values = round(gain*values) + zp;
logposrange = values > xjoint;
lognegrange = values < -xjoint;
values(logposrange) = log(values(logposrange)) + 1;
values(lognegrange) = -(log(-values(lognegrange)) + 1);
values = round(gain * values + zp);
