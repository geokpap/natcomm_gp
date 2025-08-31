function [aggval, aggpd] = dg_aggregatePasteups(vals, pds, varargin)
% For use with output from lfp_makepasteup(..., 'data').
%INPUTS
% vals: an array of <values> structures as returned by
%   lfp_makepasteup(..., 'data').
% pds: an array of <plotdata> structures as returned by
%   lfp_makepasteup(..., 'data'); must be same size as <values>.  The value
%   of  the field "aligns" must be identical in every element of <pds>.
%OUTPUTS
% aggval: a single <values> structure as returned by
%   lfp_makepasteup(..., 'data') containing the aggregation (i.e. the
%   vertical concatenation) of all the individual elements of <vals>.  It
%   is an error if different elements of <vals> contain different sets of
%   fields.
% aggpd: a single <plotdata> structure as returned by
%   lfp_makepasteup(..., 'data') with the fields "timepts", "win", and
%   "offsets" modified as needed to match <aggval>.  Fields that require no
%   modifications are copied verbatim if they are the same in all elements
%   of <pds>, or are replaced by cell vectors containing one value from one
%   element of <pds> per cell.  There is also an additional field,
%   "commonwinabs", containing the absolute time windows covered by each
%   alignment event in the final aggregated pasteup.
%OPTIONS
% 'boots' - handles *_boot cell arrays by performing the vertical
%   concatenation on the contents of each each row of the cell array.  Also
% 	repairs pasteups that contain mishandled *_boot cell arrays that have a
% 	separate column for each alignment instead of horizontally concatenated
% 	cell contents.

%$Rev: 151 $
%$Date: 2012-06-18 17:06:32 -0400 (Mon, 18 Jun 2012) $
%$Author: dgibson $

bootsflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'boots'
            bootsflag = true;
        otherwise
            error('dg_aggregatePasteups:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isequal(size(vals), size(pds))
    error('dg_aggregatePasteups:size', ...
        '<vals> and <pds> must both be the same size.');
end
if numel(vals) == 1
    aggval = vals;
    aggpd = pds;
    return
end

valfields = fieldnames(vals(1));
pdfields = setdiff(fieldnames(pds(1)), {'timepts' 'offsets'});
aggpd = pds(1);
isconstantpd = true(size(pdfields));

% Find the common time framework.  <IEIs> and <offsets> are both in
% alignment X vals format.  Also find the constant fields in <pds>.
sampleperiod = median(diff(pds(1).timepts));
aligns = pds(1).aligns;
% The entire first row of <offsets> is zero by definition, and explicitly
% assigned here for later convenience.
offsets = zeros(numel(aligns), numel(vals));
offsets(2:end, 1) = reshape(pds(1).offsets, [], 1);
IEIs = NaN(numel(pds(1).offsets), numel(vals));
IEIs(:,1) = diff(offsets(:,1));
starttimes = NaN(1,numel(vals));    % pasteup start times
starttimes(1) = pds(1).timepts(1);
endtimes = NaN(1,numel(vals));      % pasteup end times
endtimes(1) = pds(1).timepts(end);
for idx = 2:numel(vals)
    if ~isequal(pds(idx).aligns, aligns)
        error('dg_aggregatePasteups:aligns', ...
            '"aligns" must be identical in every element of <pds>');
    end
    if abs(median(diff(pds(idx).timepts)) - sampleperiod) > 1e-6
        error('dg_aggregatePasteups:period', ...
            '"timepts" must have the same median sample period in every element of <pds>');
    end
    offsets(2:end, idx) = reshape(pds(idx).offsets, [], 1);
    IEIs(:,idx) = diff(offsets(:,idx));
    starttimes(idx) = pds(idx).timepts(1);
    endtimes(idx) = pds(idx).timepts(end);
    % Check for constant fields in <pds> while we're at it:
    for fieldnum = 1:length(pdfields)
        if isconstantpd(fieldnum) && ~isequal( ...
                pds(idx).(pdfields{fieldnum}), ...
                aggpd.(pdfields{fieldnum}) )
            isconstantpd(fieldnum) = false;
        end
    end
end
% <commonwin> is in aligns X [start end] format, always
% relative to the corresponding alignment event.
commonwin = NaN(length(aligns), 2);
commonwin(1, 1) = max(starttimes);
commonwin(1:end-1, 2) = min(IEIs, [], 2) / 2;
commonwin(2:end, 1) = -commonwin(1:end-1, 2);
commonwin(end, 2) = min(endtimes - offsets(end,:));

% Now we construct an array of column indices that map from input arrays to
% output arrays.  We begin with a cell array because there will sometimes
% be mismatches in numbers of columns, which must be reconciled before we
% can convert to a plain numeric array.  In the spirit of "common time",
% reconciliation is done by deleting samples.  In order to prevent
% headaches, each element of <idxmap> is a row vector, so that columns
% represent - yes! - columns.  For easy conversion once it's reconciled,
% <idxmap> is in vals X alignment format, which is the transpose of the
% format for <IEIs> and <offsets>.  <relTS> is the timestamps of the
% selected columns relative to their alignment reference.
idxmap = cell(numel(vals), numel(aligns));
relTS = cell(numel(vals), numel(aligns));
for idx = 1:numel(vals)
    for alignidx = 1:numel(aligns)
        % Convert <commonwin> to absolute times for <commonwinabs>:
        commonwinabs = commonwin(alignidx,:) + offsets(alignidx, idx);
        % The convention is to assign a column to an event window if the
        % column's timept is within the window.  In case of equality, the
        % column goes with the window to the left, except for the left hand
        % edge of the very first window (which is included in that window). 
        if alignidx == 1
            isinwindow = pds(idx).timepts >= commonwinabs(1) & ...
                pds(idx).timepts <= commonwinabs(2);
        else
            isinwindow = pds(idx).timepts > commonwinabs(1) & ...
                pds(idx).timepts <= commonwinabs(2);
        end
        idxmap{idx, alignidx} = reshape(find(isinwindow), 1, []);
        relTS{idx, alignidx} = reshape( pds(idx).timepts( ...
            idxmap{idx, alignidx} ) - offsets(alignidx, idx), 1, [] );
        if idx > 1 && length(idxmap{idx, alignidx}) ~= ...
                length(idxmap{idx-1, alignidx})
            % Must trim off the longer one(s), but on which side?  We align
            % the new relTS to the median of the existing relTS by
            % minimizing the sum of the absolute differences.
            lengthdiff = abs( length(idxmap{idx, alignidx}) - ...
                length(idxmap{idx-1, alignidx}) );
            medianTS = median(cell2mat(relTS(1:idx-1, alignidx)), 1);
            if length(idxmap{idx, alignidx}) > ...
                    length(medianTS)
                longerTS = idxmap{idx, alignidx};
                shorterTS = medianTS;
            else
                longerTS = medianTS;
                shorterTS = idxmap{idx, alignidx};
            end
            sumdiff = NaN(lengthdiff+1,1);
            for startpt = 1:(lengthdiff+1)
                sumdiff(startpt) = sum(abs( shorterTS - ...
                    longerTS(startpt : (length(shorterTS)+startpt-1)) ));
            end
            [m, startpt] = min(sumdiff); %#ok<ASGLU>
            cols2delL = 1 : (startpt-1);
            cols2delR = ...
                (length(longerTS)-lengthdiff+startpt) : length(longerTS);
            if length(idxmap{idx, alignidx}) > ...
                    length(medianTS)
                % Delete from new entries
                idxmap{idx, alignidx}([cols2delL cols2delR]) = [];
                relTS{idx, alignidx}([cols2delL cols2delR]) = [];
            else
                % Delete from each of the old entries
                for oldidx = 1:idx-1
                    idxmap{oldidx, alignidx}([cols2delL cols2delR]) = [];
                    relTS{oldidx, alignidx}([cols2delL cols2delR]) = [];
                end
            end
        end
    end
end
% <colidx> is in vals X column-indices format:
colidx = cell2mat(idxmap);

% Now we do the vertical concatenation on every field in <vals>.
for fieldnum = 1:length(valfields)
    aggval.(valfields{fieldnum}) = [];
end
for idx = 1:numel(vals)
    for fieldnum = 1:length(valfields)
        if bootsflag && ~isempty( ...
                regexp(valfields{fieldnum}, '_boot$', 'once') )
            for bootnum = 1:size(vals(idx).(valfields{fieldnum}), 1)
                pasteupval = ...
                    cell2mat(vals(idx).(valfields{fieldnum})(bootnum,:));
                if idx == 1 
                    if bootnum == 1
                        % The first time through,
                        % aggval.(valfields{fieldnum}) will be the empty
                        % numeric [], and must be replaced by a single
                        % element cell.
                        aggval.(valfields{fieldnum}) = ...
                            {pasteupval(:, colidx(idx,:))};
                    else
                        % We are adding another boot from the first element
                        % of <vals>, so no concatenation this time.
                        aggval.(valfields{fieldnum}){bootnum} = ...
                            pasteupval(:, colidx(idx,:));
                    end
                else
                    aggval.(valfields{fieldnum}){bootnum} = [
                        aggval.(valfields{fieldnum}){bootnum}
                        pasteupval(:, colidx(idx,:)) ];
                end
            end
        else
            aggval.(valfields{fieldnum}) = [ aggval.(valfields{fieldnum})
                vals(idx).(valfields{fieldnum})(:, colidx(idx,:)) ];
        end
    end
    % ...and also handle the non-constant fields in <pds> while we're at it
    if any(~isconstantpd)
        for fieldnum = reshape(find(~isconstantpd), 1, [])
            if idx==1
                aggpd.(pdfields{fieldnum}) = {aggpd.(pdfields{fieldnum})};
            else
                aggpd.(pdfields{fieldnum}){end+1} = ...
                    pds(idx).(pdfields{fieldnum});
            end
        end
    end
end

% And finally fudge the adjusted output fields.
% The entire first row of <offsets> is zero by definition, and is thus
% omitted by convention of lfp_makepasteup.
aggpd.offsets = cumsum([0; min(IEIs, [], 2)]);
aggpd.commonwinabs = commonwin + repmat(aggpd.offsets, 1, 2);
aggpd.win = [aggpd.commonwinabs(1,1) aggpd.commonwinabs(end,2)];
absTS = cell(1, numel(aligns));
for alignidx = 1:numel(aligns)
    absTS{alignidx} = median( cell2mat(relTS(:,alignidx)) ...
        + repmat(offsets(alignidx,:)', 1, size(relTS{1,alignidx}, 2)), 1 );
end
aggpd.timepts = cell2mat(absTS);
aggpd.offsets = aggpd.offsets(2:end);
