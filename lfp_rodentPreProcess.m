function events = lfp_rodentPreProcess(events)
% Discussion with Terra and Yasuo relating to reading VT-based event 23s
% to fix certain Jianbin event files that lack them resulted in the
% following new strategy for merging events from Nlx files with events from
% VT output files.  If there is a VT events file, then all events from it
% are read and merged with the events from the Nlx file, subject to the
% constraint that no duplicate entries (i.e. with both the same event ID
% and the same timestamp) will result. The "no duplicate entries" rule is
% implemented using dg_tolUnion, so no matter how ridiculous the truncation
% error problems get, it doesn't matter.  This strategy "should" work
% equally well for both mouse and rat data.
%   - DG 25-Feb-2008

%$Rev: 280 $
%$Date: 2012-07-17 18:42:05 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;
% Rodent events use only the lower byte.  However, there are sometimes
% outboard TTL-controlled devices (laser, reward pumps) that are controlled
% by bits in the high order byte.  These might constitute events unto
% themselves (e.g. reward pumps) or they might piggyback on other events
% (laser, which can be superimposed on a series of consecutive events).  In
% the case where you want to unconditionally ignore certain bits, they can
% be masked out be setting a value for EventBitMask in your getEvtIDs file.
if ~exist('EventBitMask', 'var')
    EventBitMask = hex2dec('FFFF');
end
events(:,2) = bitand(events(:,2), EventBitMask);
[parentdir,sessionID] = fileparts(lfp_DataDir);
sessionID = upper(sessionID);
% Terra's "New Events" file:
fnameEV=[parentdir '\NE' sessionID '.dat'];
if exist(fnameEV, 'file')
    fprintf(1, 'Opening MVT events file "%s".\n', fnameEV);
    [TimeStamps, TTLIDs] = lfp_readYasuoEvents(fnameEV);
else
    % Yasuo's Delphi output file:
    fnameEV=[parentdir '\E' sessionID '.dat'];
    if exist(fnameEV, 'file')
        fprintf(1, 'Opening VT events file "%s".\n', fnameEV);
        [TimeStamps, TTLIDs] = lfp_readYasuoEvents(fnameEV);
    else
        warning('lfp_rodentPreProcess:missingVTevts', ...
            'Could not find VT events file "%s"; proceeding without.', ...
            fnameEV);
        return;
    end
end
% E<sessionID>.dat has 10kHz timestamps:
TimeStamps = 100 * TimeStamps;
% Idiot check:
if ~(TimeStamps(1)>0.99*events(1,1) && TimeStamps(end)<1.01*events(end,1))
    error('lfp_rodentPreProcess:badVTfile', ...
        'The timestamps in %s are not in the range spanned by %s!', ...
        fnameEV, lfp_DataDir );
end
events = dg_tolUnion(events(:, [2 1]), [TTLIDs' TimeStamps'], 100);
events = sortrows(events(:, [2 1]));

