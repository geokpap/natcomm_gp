function [offset, targidx, matchidx, time_est] = dg_findTrialOffset( ...
    targTS, targID, rawTS, rawIDs, rawrange, force);
% Finds an integer offset to add to the event times given in <targTS>
% so that a best match against the times in <rawTS> and events in <rawIDs>
% is obtained.  Matching is restricted to rawTS(rawrange), where <rawrange>
% is a numeric (not logical) index array.  Elements are considered to match
% if the difference between their timestamp values is < 2 and their IDs are
% equal.  If multiple elements of <rawTS> match a single element of
% <targTS>, then the closest value is used. If multiple offsets are
% tied, then multiple values are returned. <matchidx> is a cell row vector
% each of whose elements is a column vector of indices into <rawTS> of the
% elements that matched <targTS> at the corresponding offset. <targidx>
% is a corresponding list of indices into <targTS>, and is of the same
% structure as <matchidx>. Only non-negative offsets are considered.
%
% Before starting calculation, an estimate of computing time is made.  If
% <force> is [], then a dialog box is presented when the computing time is
% expected to be large.  If <force> = 1 then the computation proceeds
% regardless, and if <force> = 0 then the computation is abandoned when
% expected computing time is large.  In the latter case, empty values are
% returned for everything except <time_est>, which contains the
% time estimate in minutes, or the value -1 if execution was aborted for
% some other reason.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

offset = [];
targidx = cell(1,0);  
matchidx = cell(1,0);  
time_est = -1;
if isempty(rawrange)
    return
end
rawTS2 = reshape(rawTS(rawrange), [], 1);
rawIDs2 = reshape(rawIDs(rawrange), [], 1);
maxoffset = round(max(rawTS2) - min(targTS));
minoffset = max(0, round(min(rawTS2) - max(targTS)));
%minutesPerComparison = 1.4e-009; % Dell Precision M65 1.66 GHz
a = 8.3e-008;  % minutes, Dell Precision M65 1.66 GHz
b = 1.1e-009;  % minutes, Dell Precision M65 1.66 GHz
time_est = (maxoffset-minoffset)*length(targTS) * (a + b*length(rawrange));
disp(sprintf('%d:%d time_est=%d:', rawrange(1), rawrange(end), time_est));
disp(sprintf( ...
    'maxoffset=%d; minoffset=%d; #evts=%d; length(rawrange)=%d',...
    maxoffset, minoffset, length(targTS), length(rawrange) ));
if time_est > 20
    if isempty(force)
        str = sprintf( ...
            'This computation will take %d minutes at 1.66 GHz; continue?', ...
            time_est );
        button = questdlg(str,'Run long computation?','Yes','No','Details','Yes');
        switch button
            case 'No'
            return
            case 'Details'
                str = sprintf( ...
                    'max(rawTS2): %d min(rawTS2): %d\nmax(targTS): %d min(targTS): %d', ...
                    max(rawTS2), min(rawTS2), max(targTS), min(targTS) );
                button = questdlg(str,'Run long computation?','Yes','No','Yes');
                if isequal(button, 'No')
                    return
                end
        end
    else
        switch force
            case 0
                return
            case 1
                % do nothing, just keep going
            case {2 3}
                % for 2006-3-17_16-0-12	nickel 0301706.1, use first 100
                % events.
                if force==2
                    rawrange = rawrange(1) : ...
                        min(rawrange(1)+99,length(rawrange));
                else % case == 3
                    rawrange = max(1,length(rawrange)-99):rawrange(end);
                end
                if isempty(rawrange)
                    disp('Recalculated rawrange is empty');
                    return
                else
                    rawTS2 = reshape(rawTS(rawrange), [], 1);
                    rawIDs2 = reshape(rawIDs(rawrange), [], 1);
                    maxoffset = round(max(rawTS2) - min(targTS));
                    minoffset = max(0, round(min(rawTS2) - max(targTS)));
                    time_est = (maxoffset-minoffset)*length(targTS) * ...
                        (a + b*length(rawrange));
                    disp(sprintf('Invoked force=%d; %d:%d time_est=%d:', ...
                        force, rawrange(1), rawrange(end), time_est));
                    disp(sprintf( ...
                        'maxoffset=%d; minoffset=%d; #evts=%d; length(rawrange)=%d',...
                        maxoffset, minoffset, length(targTS), length(rawrange) ));
                end
        end
    end
end
for shift = minoffset:maxoffset
    idxlist = zeros(0,1);
    trialevtidxlist = zeros(0,1);
    for trialevtidx = 1:length(targTS)
        idx = find( ...
            (abs(rawTS2 - targTS(trialevtidx) - shift) < 2) ...
            & (rawIDs2 == targID(trialevtidx)) );
        if length(idx) > 1
            [c, i] = min( ...
                abs(rawTS2(idx) - targID(trialevtidx)  - shift) );
            idx = idx(i);
        end
        if ~isempty(idx)
            idxlist(end+1) = idx;
            trialevtidxlist(end+1) = trialevtidx;
        end
    end
    if isempty(matchidx) || length(idxlist) > length(matchidx{1})
        % This is a better match
        matchidx = {idxlist};
        targidx = {trialevtidxlist};
        offset = shift;
    elseif length(idxlist) == length(matchidx{1})
        % This is a tie
        matchidx(end+1) = {idxlist};
        targidx(end+1) = {trialevtidxlist};
        offset(end+1) = shift;
    end
end
for k = 1:length(matchidx)
    if ~isequal(size(targidx{k}), size(matchidx{k}))
        error('dg_findTrialOffset:idxsync', ...
            'Internal error: targidx and matchidx are out of sync.' );
    end
    % matchidx points into rawTS2, but we need to return values that point
    % into rawTS:
    matchidx{k} = matchidx{k} + rawrange(1) - 1;
end
% Eliminate degenerate solutions, defined as solutions where targidx and
% matchidx are the same.  First, find them. degenidxlist will contain each
% member of each set of degenerate solutions, and degensets will contain
% the members of each set:
degenidxlist = zeros(1,0);
degensets = cell(1,0);
for k = 1:length(offset)
    if ismember(k, degenidxlist)
        continue
    end
    degenidx = k;
    for j = (k+1):length(offset)
        if isequal(matchidx{j}, matchidx{k}) ...
                && isequal(targidx{j}, targidx{k})
            degenidx = [degenidx j];
        end
    end
    if length(degenidx) > 1
        degenidxlist = [degenidxlist degenidx];
        degensets(end+1) = {degenidx};
    end
end
% Then, for each set of degenerate solutions, remove all but the one that
% yields the smallest absolute timestamp error:
deleteidx = [];
for k = 1:length(degensets)
    errors = zeros(1,0);
    for degenidx = degensets{k}
        errors(end+1) = mean(abs( ...
            reshape(rawTS(matchidx{degenidx}), 1, []) ...
            - reshape(targTS(targidx{degenidx}), 1, []) ));
    end
    [c,i] = min(errors);
    deleteidx = [deleteidx setdiff(degensets{k}, degensets{k}(i))];
end
offset(deleteidx) = [];
targidx(deleteidx) = [];
matchidx(deleteidx) = [];
       
    
