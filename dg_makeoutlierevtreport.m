function [report] = dg_makeoutlierevtreport( ...
    eventlist, allintervals, tailprct, alltrialIDs, detailflag)

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

report = sprintf('Outlier Interevent Intervals\n');
tails = prctile(allintervals, [tailprct 100-tailprct]);
for evtidx = 1:(length(eventlist)-1)
    outliers = allintervals(:,evtidx) < tails(1, evtidx) | ...
        allintervals(:,evtidx) > tails(2, evtidx);
    intervalstring = sprintf('%s - %s', ...
        mat2str(eventlist{evtidx}), mat2str(eventlist{evtidx+1}) );
    report = sprintf('%s%-12s:\t\n', ...
        report, intervalstring );
    if detailflag
        for outlieridx = find(outliers')
            report = sprintf('%s\t\t%s\n', report, alltrialIDs{outlieridx});
        end
    end
end