function goodttl = dg_plotTMazeEvtSeq(filename, varargin)
%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

jianbinflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'jianbin'
            jianbinflag = true;
        otherwise
            error('funcname:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% evtorder<0 is for housekeeping events such as rec on/off, trial outcome,
% etc. isnan(evtorder) implies that we aren't even looking at it.
evtorder = NaN(32767, 1);
evtorder(1) = -1;
evtorder(2) = -2;
evtorder(90:99) = -3;
evtorder(11) = 2;   % Gate
evtorder(13) = 3;   % Out of Start
evtorder(6) = 4;    % Pre-tone
evtorder(9) = 6;    % Pre-turn
evtorder(14) = 7;   % Turn
evtorder(15) = 8;   % Turn Off
evtorder(16) = 8;   % Turn Off
evtorder(7) = 9;    % Pre-goal
evtorder(8) = 9;    % Pre-goal
evtorder(17) = 10;  % Goal
evtorder(18) = 10;  % Goal
evtorder(41) = 11;  % Tone Off
evtorder(48) = 11;  % Tone Off
if jianbinflag
    evtorder(31) = 1;   % Tone On
    evtorder(38) = 1;   % Tone On
    evtorder(23) = 5;   % Mid-stem ("Tone On" photobeam)
else
    evtorder(10) = 1;   % Click
    evtorder(31) = 5;   % Tone On
    evtorder(38) = 5;   % Tone On
    evtorder(21) = 5;   % Tactile Cue
    evtorder(22) = 5;   % Tactile Cue
end

[TS, TTL, ES, Hdr] = dg_readEvents(filename);
isgoodstrobe = TTL<0;
goodttl = reshape(TTL(isgoodstrobe) + 2^15, [], 1);
ttlorder = evtorder(goodttl);
hF = figure;
plot(ttlorder, '-+');
hold on;
plot(get(gca, 'XLim'), [0 0], 'k');
ylim([-20 20]);



