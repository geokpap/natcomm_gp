function evtfname = dg_findEventsFile(sessiondir, evtExprs)
% evtfname = dg_findEventsFile(sessiondir)
% evtfname = dg_findEventsFile(sessiondir, evtExprs)
%INPUTS
% sessiondir: directory to search for an events file.
% evtExprs: this argument is optional, and specifies a list of regular
%   expressions to use in order of decreasing preference to identify an
%   events file.  See code for default.

%$Rev: 239 $
%$Date: 2016-03-10 19:30:16 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $

if nargin < 2
    evtExprs = {'events.*\.evtsav', 'events\.nev', 'events\.mat', ...
        'events\.dat'};
end
files = dir(sessiondir);
filenames = reshape({files(~cell2mat({files.isdir})).name}, [], 1);
for evtXidx = 1:length(evtExprs)
    evtmatchExpr = evtExprs{evtXidx};
    fnidx = find(~cellfun(@isempty, regexpi(filenames, evtmatchExpr)), 1);
    if ~isempty(fnidx)
        break
    end
    if evtXidx == length(evtExprs)
        error('lfp_markBadCSC:noevents', ...
            'No events file in %s', sessiondir);
    end
end
evtfname = filenames{fnidx};

