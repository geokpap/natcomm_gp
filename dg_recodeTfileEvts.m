function rc = dg_recodeTfileEvts(infile, outfile, evtmap, varargin)
%   Reads the T file at <infile>, changes the event IDs as specified by
%   <evtmap>, and writes the result to <outfile>.  The first column of
%   <evtmap> specifies the input event IDs, and the second column specifies
%   the corresponding output event IDs.  Timestamps are moved from the
%   input ID to the output ID, i.e. the timestamps of the input IDs are set
%   to zero in the output file (unless they get overwritten by other moved
%   timestamps).  Event IDs that are not listed at all in <evtmap> are left
%   unaltered.  If <outfile> already exists, it is silently overwritten
%   (unless it is read-only, in which case the write fails and raises an
%   error).  If the input file format is 0xFFF1, then the output is written
%   in format 0xFFF3, and the missing fields are filled in by formula.
%   <rc> is true if a new file was written (see OPTIONS), and false
%   otherwise.
%OPTIONS
%   'HasEvent', evtIDs - <evtIDs> is a vector of numeric event IDs.  The T
%       file must contain at least one of the events in <evtIDs> in order
%       to be recoded; otherwise a warning is raised and nothing else
%       happens.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

rc = false;
if ~isa(evtmap, 'numeric')
    error('dg_recodeTfileEvts:evtmap1', ...
        '<evtamp> must be numeric' );
end
if isempty(evtmap)
    error('dg_recodeTfileEvts:evtmap2', ...
        '<evtamp> is empty.  Use copyfile if you want to copy the file.' );
end

evtIDs = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'HasEvent'
            argnum = argnum + 1;
            evtIDs = varargin{argnum};
        otherwise
            error('dg_recodeTfileEvts:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[fileheader, trialdata] = ...
    dg_ReadRodentFormat(infile);
if fileheader.TSize == 0
    error('dg_recodeTfileEvts:notrials', ...
        'The T file %s contains zero trials.', infile);
end
eventses = cell(size(trialdata));
[eventses{:}] = deal(trialdata.events);
eventses = cell2mat(eventses'); % eventses is in trial X evtID format
if ~isempty(evtIDs)
    if ~(any(any(eventses(:,evtIDs))))
        warning('dg_recodeTfileEvts:noevtIDs', ...
            'The file %s does not contain any of the events %s.', ...
            infile, dg_thing2str(evtIDs) );
        return
    end
end
neweventses = eventses;
neweventses(:, evtmap(:,1)) = 0;
neweventses(:, evtmap(:,2)) = eventses(:, evtmap(:,1));
neweventses = mat2cell(neweventses, ...
    ones(1,size(neweventses,1)), size(neweventses,2) );
[trialdata.events] = deal(neweventses{:});
if fileheader.Format == 65521
    fileheader.Format = 65523;
    fileheader.ProcType = 1;
    fileheader.NoGoCS = 99;
    if fileheader.RightCS == 1
        fileheader.LeftCS = 8;
    else
        fileheader.LeftCS = 1;
    end
    for trial = 1:fileheader.TSize
        trialdata(trial).SSize(11:12) = 0;
        trialdata(trial).header.SSize(11:12) = 0;
    end
end
dg_writeRodentFormat(outfile, fileheader, trialdata);
rc = true;


