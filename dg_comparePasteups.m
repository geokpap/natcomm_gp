function dg_comparePasteups(hF1, hF2, varargin)
% Modifies the line plot figure whose handle is <hF1> to show confidence
% interval comparison with <hF2>.  Either figure could be output from
% lfp_makepasteup or directly from lfp_disp.  Specifically:
% (1) copy the event markers from hF2 to hF1, changing the line style to
% dashed
% (2) change the color of the confidence bounds on the hF1 plot to red
% (3) copy the confidence bounds from the hF2 plot to the hF1 plot,
% changing their color to blue
% Note that the figure <hF2> may also be modified when using certain
% options, e.g. 'nodetails'.
% (4) append the title of hF2 to that of hF1 with "vs" in between
% (5) append any 'UserData' from <hF2> to that in <hF1>.
%
%OPTIONS
% 'alignevents' - instead of copying event markers from hF2 to hF1,
%   align the event markers and line plots as follows.  First, find the
%   time window occupied by each line segment in hF1 and hF2, and for each
%   line segment, compute the time window that is common to hF1 and hF2.
%   Horizontally shift the event markers in hF1 as needed to display only
%   the common time window around each marker.  Truncate and horizontally
%   shift each line segment in hF1 and each one copied from hF2 to match
%   the adjusted event markers.  Any common time window that would be fewer
%   than three points wide is eliminated and raises a warning.  Throws an
%   error if there are different numbers of line segments in hF1 and hF2,
%   different numbers of event markers in hF1 and hF2, or different numbers
%   of line segments and event markers.
% 'baselines', [b1 b2] - subtracts <b1> from all waveform vertical
%   coordinates in hF1 and <b2> from those in hF2.  Note that this does not
%   affect the vertical positions of any other graphics objects such as
%   even
% 'yscale', [s1 s2] - multiplies all waveform vertical
%   coordinates in hF1 by s1 and those in hF2 by s2. If used with
%   'baseslines' option the baseline-subtracted waveforms are scaled.
% 'color1', colorspec1 - overrides the default red color of waveform
%   segments in hF1; <colorspec> can be a character or a color
%   triple.
% 'color2', colorspec2 - overrides the default blue color of waveform
%   segments copied from hF2; <colorspec> can be a character or a color
%   triple.
% 'linewidth', linewidth - sets wave plots line width to <linewidth>.  For
%   no good reason, this only works in combination with 'alignevents'.
% 'nodetails' - uses only those events whose value for 'ButtonDownFcn' is
%   'fprintf(1,'No details available\n')'.  If the following event is of
%   the same color and within 1 ms of the 'nodetails' event, then its
%   'ButtonDownFcn' is copied.  If not, a warning is raised.
% 'rmavg' - delete the mean data line from the hF1 plot
% 'nocopy', ts - <ts> is a vector of time coordinates in seconds; skips
%   copying any event that occurs within 1 millisecond of any member of
%   <ts>.
% 'dashedevents' - if 'alignevents' option is used, makes all event markers dashed,
% in order to aid visualization. Only works with 'alignevents'.

%$Rev: 243 $
%$Date: 2016-06-24 21:51:26 -0400 (Fri, 24 Jun 2016) $
%$Author: dgibson $

lfp_declareGlobals;

aligneventsflag = false;
linewidth = [];
nodetailsflag = false;
rmavgflag = false;
ts = [];
argnum = 1;
baselines = [];
scalefactor = [];
colorspec1 = 'r';
colorspec2 = 'b';
dashedeventsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alignevents'
            aligneventsflag = true;
        case 'baselines'
            argnum = argnum + 1;
            baselines = varargin{argnum};
        case 'yscale'
            argnum = argnum + 1;
            scalefactor = varargin{argnum};
        case 'color1'
            argnum = argnum + 1;
            colorspec1 = varargin{argnum};
        case 'color2'
            argnum = argnum + 1;
            colorspec2 = varargin{argnum};
        case 'linewidth'
            argnum = argnum + 1;
            linewidth = varargin{argnum};
        case 'nocopy'
            argnum = argnum + 1;
            ts = varargin{argnum};
        case 'nodetails'
            nodetailsflag = true;
        case 'rmavg'
            rmavgflag = true;
        case 'dashedevents'
            dashedeventsflag = true;
        otherwise
            error('dg_comparePasteups:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Get hF1 data
[evts1, waves1, extraevts] = dg_getlinedata(hF1, nodetailsflag);
for hL = extraevts'
    delete(hL);
end

% if 'alignevents' option, make event markers of hF1 dashed, so they interfere
% less with the visualization of the data curves:
if aligneventsflag && dashedeventsflag
    for hL = evts1
        if hL
            set(hL, 'LineStyle', '--');
        end
    end
end

% Change hF1 coloring
for hL = waves1(1:2,:)
    if hL
        set(hL, 'Color', colorspec1);
    end
end

hA1 = get(hF1, 'CurrentAxes');
ylimits = get(hA1, 'YLim');

% Get hF2 data
[evts2, waves2] = dg_getlinedata(hF2, nodetailsflag);
hA2 = get(hF2, 'CurrentAxes');

if aligneventsflag
    % The line segments are in reverse order when the plots are produced by
    % lfp_makePasteup, so if that is the case, reverse the order of the
    % wave segments and events markers lists.
    if isreverse(evts1)
        evts1 = evts1(end:-1:1);
        waves1 = waves1(:,end:-1:1);
    end
    if isreverse(evts2)
        evts2 = evts2(end:-1:1);
        waves2 = waves2(:,end:-1:1);
    end
    ct = findcommontime(evts1, waves1, evts2, waves2);
    % shift event markers
    for k = 1:length(evts1)
        adjustlinesegment(ct(1), evts1(k), []);
    end
    % truncate and shift hF1 line segments
    for rownum = 1:size(waves1,1)
        for hL = waves1(rownum,:)
            if hL
                adjustlinesegment(ct(1), hL, linewidth);
            end
            if ~isempty(baselines)
                set(hL, 'YData', get(hL, 'YData') - baselines(1));
            end
            if ~isempty(scalefactor)
                set(hL, 'YData', get(hL, 'YData') * scalefactor(1));
            end
        end
    end
else
    % Copy event markers
    for hL = evts2'
        xdata = get(hL, 'XData');
        if isempty(ts) || all(abs(ts - xdata(1)) > 1e-3)
            hCopy = copyobj(hL, hA1);
            set(hCopy, 'LineStyle', '--');
            set(hCopy, 'YData', ylimits);
        end
    end
end

% Copy line plot segments
for k = 1:2
    for hL = waves2(k,:)
        if hL
            hCopy = copyobj(hL, hA1);
            set(hCopy, 'Color', colorspec2);
            if aligneventsflag
                adjustlinesegment(ct(2), hCopy, linewidth);
            end
            if ~isempty(baselines)
                set(hCopy, 'YData', get(hCopy, 'YData') - baselines(2));
            end
            if ~isempty(scalefactor)
                set(hCopy, 'YData', get(hCopy, 'YData') * scalefactor(2));
            end
        end
    end
end

% Remove avg curves if requested
if rmavgflag
    for hL = waves1(3,:)
        if hL
            delete(hL);
        end
    end
end

% Append titles
str1 = get(get(hA1, 'Title'), 'String');
if iscell(str1)
    tmpstr = str1{1};
    for k = 2:length(str1)
        tmpstr = sprintf('%s\n%s', tmpstr, str1{k});
    end
    str1 = tmpstr;
end
str2 = get(get(hA2, 'Title'), 'String');
if iscell(str2)
    tmpstr = str2{1};
    for k = 2:length(str2)
        tmpstr = sprintf('%s\n%s', tmpstr, str2{k});
    end
    str2 = tmpstr;
end
hT = title(hA1, sprintf('%s vs %s\n', str1(1,:), str2(1,:)), 'Interpreter', 'none');

% Check hF2 for UserData
userdata2 = get(hF2, 'UserData');
if ~isempty(userdata2)
    userdata1 = get(hF1, 'UserData');
    if isempty(userdata1)
        set(hF1, 'UserData', userdata2);
    else
        set(hF1, 'UserData', ...
            [reshape(userdata1, 1, []) reshape(userdata2, 1, [])] );
    end
    set(hT, 'ButtonDownFcn', ...
        'lfp_makepasteup_userdata=get(gcf,''UserData''); disp(''UserData in lfp_makepasteup_userdata''); disp(dg_thing2str(lfp_makepasteup_userdata));' );
end
end


function ct = findcommontime(evts1, waves1, evts2, waves2)
% evts1, waves1, evts2, waves2 must all contain the same number of objects
% (one per alignment event aka time window).
% <ct> contains all the data needed to truncate and shift line objects
% whether they were originally in hF1 or hF2.  ct(1) refers to objects
% originating in hF1, ct(2) to hF2.  The fields all have one row for
% each line segment, and are:
%   intervals:  a 2-column array of start and end times
%   truncleft:  column vector of numbers of points to remove from the
%       beginning
%   truncright:  numbers of points to remove from the end
%   shift:  values to add to old x-coordinates to shift them to new
%       positions
% There are two incremental contributions to each shift:  the time deleted
% from the right side of the preceding segment, and the time deleted from
% the left side of the current segment.  The total shift for any given
% segment is thus the running total of the incremental shifts summed
% through the given segment.  The first segment is by definition time zero
% and is not to be shifted.  There is also an overriding constraint on the
% shifts for the two plots, namely that after the shifts have been made,
% the event markers must be shifted to exactly the same time, so only one
% figure can have its shifts calculated by summing its incremental shifts.
if length(evts1) ~= length(evts2)
    error('dg_comparePasteups:evtmismatch', ...
        'The plots contain different numbers of event markers.' );
elseif size(waves1,2) ~= size(waves2,2)
    error('dg_comparePasteups:wavemismatch', ...
        'The plots contain different numbers of line segments.' );
elseif size(waves1,2) ~= length(evts1)
    error('dg_comparePasteups:mismatch', ...
        'The plots contain different numbers of line segments and events.' );
end
ct(1).shift = zeros(size(evts1));
ct(2).shift = zeros(size(evts1));
for k = 1:length(evts1)
    x = get(evts1(k), 'XData');
    et1(k) = x(1); %#ok<AGROW>
    x = get(evts2(k), 'XData');
    et2(k) = x(1); %#ok<AGROW>
end
for evtix = 1:length(evts1)
    ts1 = get(waves1(3,evtix), 'XData');
    T1 = median(diff(ts1));
    ct(1).intervals(evtix,:) = ts1([1 end]);
    ts2 = get(waves2(3,evtix), 'XData');
    T2 = median(diff(ts2));
    ct(2).intervals(evtix,:) = ts2([1 end]);
    % find common time window relative to the alignment event
    commonwin(1) = max(ts1(1) - et1(evtix), ts2(1) - et2(evtix));
    commonwin(2) = min(ts1(end) - et1(evtix), ts2(end) - et2(evtix));
    ct(1).truncleft(evtix,1) = sum(ts1 < et1(evtix) + commonwin(1) - T1/2);
    ct(2).truncleft(evtix,1) = sum(ts2 < et2(evtix) + commonwin(1) - T2/2);
    ct(1).truncright(evtix,1) = sum(ts1 > et1(evtix) + commonwin(2) + T1/2);
    ct(2).truncright(evtix,1) = sum(ts2 > et2(evtix) + commonwin(2) + T2/2);
    if evtix > 1
        % calculate current segment's running total shifts:
        ct(1).shift(evtix) = ct(1).shift(evtix - 1) ...
            - prevRdel1 - (ts1(ct(1).truncleft(evtix,1) + 1) - ts1(1));
        ct(2).shift(evtix) = et1(evtix) + ct(1).shift(evtix) - et2(evtix);
    end
    prevRdel1 = ct(1).truncright(evtix,1) * T1;
    prevRdel2 = ct(2).truncright(evtix,1) * T2;
end
end


function adjustlinesegment(ct, hL, linewidth)
% <ct> is a scalar of the struct returned by findcommontime.
ts = get(hL, 'XData');
foundevtix = false;
for evtix = 1:length(ct.shift)
    if ts(1) >= ct.intervals(evtix,1) && ts(end) <= ct.intervals(evtix,2)
        foundevtix = true;
        break
    end
end
if ~foundevtix
    error('dg_comparePasteups:mismatch2', ...
        'This line object does not match any interval.' );
end
set(hL, 'XData', get(hL, 'XData') + ct.shift(evtix));
if length(ts)>2
    if ~isempty(linewidth)
        set(hL, 'LineWidth', linewidth);
    end
    % hL is a wave plot; truncate as needed
    if ct.truncleft(evtix) || ct.truncright(evtix)
        x = get(hL, 'XData');
        y = get(hL, 'YData');
        commonrange = 1 + ct.truncleft(evtix) : length(x) - ct.truncright(evtix);
        if length(commonrange) < 3
            warning('dg_comparePasteups:rangetoosmall', ...
                'The common time interval %.6f - %.6f (%d points) is being deleted.', ...
                ct.intervals(evtix,1), ct.intervals(evtix,2), ...
                length(commonrange) );
            delete(hL);
        else
            x = x(commonrange);
            y = y(commonrange);
            set(hL, 'XData', x);
            set(hL, 'YData', y);
        end
    end
end
end


function result = isreverse(evts)
% Returns true if the events are in reverse chronological order.
for k = 1:length(evts)
    x = get(evts(k), 'XData');
    ts(k) = x(1); %#ok<AGROW>
end
result  = all(ts(1:end-1) >= ts(2:end));
end

