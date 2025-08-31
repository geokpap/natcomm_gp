function report = dg_formatEventsReport(missingevts, disordered, varargin)
% missingevts, disordered are as returned by lfp_eventAnalysis.
%OPTIONS
% 'nodetail' - suppresses listing of trial IDs after summary counts.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

detailflag = true;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'nodetail'
            detailflag = false;
        otherwise
            error('dg_formatEventsReport:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

report = sprintf('Missing Events\n');
report = [report cellformat(missingevts, detailflag)];
report = [report sprintf('\nDisordered Events\n')];
report = [report cellformat(disordered, detailflag)];
end

function result = cellformat(cellary, detailflag)
result = '';
for k = 1:size(cellary,1)
    result = sprintf('%s%-8s:\t%.0f\n', ...
        result, cellary{k,1}, cellary{k,2} );
    if detailflag
        for j = 1:length(cellary{k,3})
            result = sprintf('%s\t\t%s\n', result, cellary{k,3}{j});
        end
    end
end
end
