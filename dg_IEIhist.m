function [n, xout, IEI] = dg_IEIhist(events, evtid1, evtid2)
%[n, xout, IEI] = dg_IEIhist(events, evtid1, evtid2)
% For each event with <evtid1>, compute the time interval to the
% immediately following <evtid2>, and histogram the result.  <n> and <xout>
% are as for Matlab 'hist' function.  <events> contains timestamps in first
% col, event IDs in second col.  <IEI> is the vector of inter-event
% intervals that went into the histogram.

%$Rev: 40 $
%$Date: 2009-10-06 19:30:18 -0400 (Tue, 06 Oct 2009) $
%$Author: dgibson $

evtidx = find(events(:,2) == evtid1);
IEI = NaN(size(evtidx));
for k = 1:length(evtidx)
    idx1 = evtidx(k);
    idx2 = idx1 + 1;
    while idx2 <= size(events,1)
        if events(idx2, 2) == evtid2
            IEI(k) = events(idx2, 1) - events(idx1, 1);
            break
        else
            idx2 = idx2 + 1;
        end
    end
    if idx2 > size(events,1)
        warning('dg_IEIhist:noevt2', ...
            'Failed to find evtid2 following evtid1 at timestamp %.6f', ...
            events(idx1, 1) );
    end
end
[n, xout] = hist(IEI);
figure;
bar(xout, n);
