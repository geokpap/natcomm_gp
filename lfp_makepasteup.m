function varargout = lfp_makepasteup(aligns, funch, varargin)
%[hOutFig, hOutAx] = lfp_makepasteup(aligns, funch)
%[values, plotdata] = lfp_makepasteup(... 'pasteup_opts', {'data'})
%   Accepts a list of alignment events and computes the median inter-event
%   intervals for each successive pair. For each cell in the list of
%   alignment events, computes result of evaluating <funch> for that
%   alignment event, in a time window extending half the median time to the
%   successive alignment event on each side.  Uses lfp_XLimAll as the left
%   side of the first window and the right side of the last window.  The
%   time window is implemented by temporarily assigning a value to
%   lfp_XLimAll (which is restored to its original value before exit),
%   which means that if <funch> is capable of overriding the value of
%   lfp_XLimAll (as most of the lfp_lib analytic functions can), then final
%   pasteup will contain overlaid plots whose widths are determined by
%   <funch>, NOT by the inter-event spacing. Creates a pasteup plot, with
%   white line markers at the junctions between image subplots. Records all
%   input parameters to <funch> in the 'UserData' field of the finished
%   fig; note that relevant global values must be added using the 'aux'
%   option.  All the values are displayed when you click on the figure
%   title, and thereafter the cell array <lfp_makepasteup_userdata>
%   contains the arguments in its first cell and the auxdata in its second.
%   Labels the alignment event markers in the pasteup with names
%   listed in lfp_EventNames.
%
%NOTES
% 1.  This function modifies and then restores the values of
%   lfp_AlignmentRef and lfp_XLimAll as part of its normal functioning.
%   This means that if it crashes instead of exiting normally, then the
%   original values of those globals will not be restored.
% 2.  To minimize both computation time and development time, I decided to
%   make the output consist of separate sets of plot objects for each
%   alignment, all cohabitating in different time spans in one 'axes'.
%   This could potentially cause headaches for graphicists and malfunctions
%   for other scripts that operate on the output figure.  Both problems
%   could be rectified by adding some moderately compute-intensive
%   interpolation operations to pack all the different alignment windows
%   into images or plot objects that span the full width of the figure.
% 3.  I have not tested this version to see if it actually works with
%   lfp_spikeAnalysis output, but "it should" (caveat emptor).
%
%INPUTS
% aligns: a cell array of lists of alternative alignment events (same
%   format as 'aligns' in lfp_spikeAnalysis).
% funch: a function handle to a function that produces a figure and returns
%   the handle to that figure as its first return value.  The figure must
%   contain an image plot, a line plot, or a patch plot.  Image plot
%   figures may contain more than one image, in which case the plot image
%   must be the last in the list returned by findobj(gcf, 'Type', 'image')
%   (as it will be for figures produced by calling imagesc first and then
%   colorbar).  Line plot figures must not contain any images or patches.
%   Note that some functions (e.g. lfp_his) may require an option such as
%   'returnfig' to force them to create the actual figure in spite of the
%   fact that the return values are being used.
% varargin: this represents all of the arguments required by <funch> when
%   it is invoked directly; please refer to the documentation for <funch>.
%   Note that some care must be taken in selecting arguments here, since
%   some combinations (e.g. those that produce multiple plots) do not make
%   sense and will yield unpredictable, usually undesirable, results.  It
%   is assumed that the first numerical argument ("numerical" includes [])
%   to <funch> is a list of trials as used by lfp_disp and lfp_spec.  Note
%   that this makes it impossible to specify trials as Unique Trial IDs
%   (which are strings), leading to unpredictable behavior; instead, use
%   lfp_getTrialNum to convert trial IDs to trialnums.
%
%OUTPUTS
% <hOutFig> - pasteup figure handle
% <hOutAx> - pasteup axes handle
%
%OPTIONS that are passed through to <funch> but also have effects in
% lfp_makepasteup:
% 'norefOK' - independently for each alignment reference, simply skips any
%   trials that are missing the alignment event.
%
%OPTIONS for lfp_makepasteup:
% 'pasteup_opts' - this is the only option that does NOT get passed on to
% <funch>.
% IMPORTANT: All of the following options must be invoked by specifying the
% lfp_makepasteup option 'pasteup_opts' followed by a cell array
% containing the options to be invoked, e.g.
%   lfp_makepasteup({10 [31 38] 50}, @lfp_disp, ...
%       [], filenum, 'avg', 'err2', 'evtavg', [11 14], ...
%       'pasteup_opts', {'bg', {[.8 .8 .8] [1 .8 .8]}});
% 'pasteup_opts' does not have to follow the arguments to <funch>, but can
% be placed anywhere in the argument list after <funch>.  Note that this
% implies that none of the arguments or options to <funch> can have the
% value 'pasteup_opts'.
% 'autoselect' - automatically de-selects any trials that do not have all
%   of the alignment events.  IMPORTANT:  this changes the state of
%   lfp_SelectedTrials!  That may or may not be a desirable side effect.
% 'aux', auxdata - includes <auxdata> in the list of parameters saved in
%   the figure.  <auxdata> may be any Matlab value, including cell and
%   string.
% 'bg', bg - <bg> is a cell array whose first element is an RGB triple
%   specifying the color of background rectangles for the odd-numbered
%   subplots in a line graph pasteup, and whose second element does
%   likewise for even-numbered plots.  If there is no second element, then
%   no background rectangle is drawn for even-numbered subplots.
% 'blankevtlabels' - instead of creating 'evt%d' labels for event IDs that
%   have empty entries in lfp_EventNames, leaves those event markers
%   unlabeled.
% 'data' - operates on and returns data structures instead of figures.  In
%   this case, <funch> must be a handle to a function whose first two
%   return values are (1) a structure (designated here as <values>)
%   containing various values whose columns represent time points and
%   either are row vectors or 2-D arrays; (2) a structure (here designated
%   <plotdata>) containing a field named "timepts" which contains a row
%   vector of time values corresponding to the columns of the fields in
%   <values>.  <values> may not contain any other fields.  <plotdata> may
%   contain additional fields, but it is assumed that all of these fields
%   will either be the same regardless of the alignment event, or will not
%   be used in constructing the final plot.  In any case, only the values
%   generated by the call to <funch> for the first alignment event are
%   passed through to the output <plotdata>.  The value of
%   <plotdata.timepts> is the concatenation of the values for each
%   alignment event.  Note that at the time points where the data from
%   different alignments are spliced together, this may result in
%   successive time points being separated by a shorter or longer interval
%   than the others.  <plotdata> as returned will have three additional
%   fields (or any pre-existing fields of the same names will be replaced),
%   namely "win", "offsets" and "aligns".  "aligns" is exactly a copy of
%   the input argument <aligns>.  "offsets" contains the times of each
%   alignment event except for the first, as offset in the pasteup,
%   expressed relative to the first alignment event.  "win" is the value of
%   lfp_XLimAll at the time of invocation of lfp_makepasteup. <values> as
%   returned in the output contains all the same fields as <values>
%   returned by <funch>, but their contents is the concatenation of the
%   contents returned for each alignment event.
% 'intervals', IEIs - overrides the calculation of median inter-event
%   intervals; <IEIs> must be a vector of inter-event intervals in seconds,
%   with one fewer elements than the number of elements in <aligns>.
% 'joinmarkers' - shows white line markers at the junctions between image
%   subplots
% 'nodisplay' - Specifying 'nodisplay' still creates the figure, but with
%   'Visible' set to 'off'.  This is convenient for some bulk processing
%   applications.  'nodisplay' must be grouped next to (either before or
%   after) 'pasteup_opts', and it will be included in 'pasteup_opts'
%   regardless.  Note that the figure continues to exist, take up storage,
%   take up its figure number, etc., even though there is no visible
%   manifestation of the figure on the screen.  Such figures should at some
%   point be closed by doing close(hOutFig).  Invisible figures that may
%   have gotten lost along the way can be found by doing
%       findobj(0, 'Type', 'figure')
%   and the numbers that are returned by the findobj call can be used as
%   arguments to "close", or anything else that requires a figure handle,
%   e.g.: "close(37)" or "set(37,'Visible','on')".  Invisible figures will
%   also be closed by doing "close all".  Note that if you save an
%   invisible figure to a file, then when you open the file the figure will
%   still be invisible, and you will have to figure out its figure number
%   and set its 'Visible' property to 'on' to see it.
% 'title', title - substitutes <title> for the default figure title.

%$Rev: 425 $
%$Date: 2023-10-20 17:24:30 -0400 (Fri, 20 Oct 2023) $
%$Author: dgibson $

lfp_declareGlobals;

if ~iscell(aligns)
    error('lfp_makepasteup:aligns', ...
        '<aligns> must be a cell array.');
end
if ~isa(funch, 'function_handle')
    error('lfp_makepasteup:funch', ...
        '<funch> must be a function handle.');
end

auxdata = {};
bg = {};
blankevtlabels = false;
norefOKflag = false;
showjoinmarkers = false;
titlestr = func2str(funch);

trials = [];
trialargnum = 0;
for k = 1:length(varargin)
    if isnumeric(varargin{k})
        trials = varargin{k};
        trialargnum = k;
        break
    end
end

% Find and edit out the pasteup_opts (if any), and handle options that
% affect both lfp_makepasteup and <funch>:
args2delete = [];
pasteup_opts = [];
argnum = 1;
display_str = 'on';
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'norefOK'
                norefOKflag = true;
            case 'pasteup_opts'
                argnum = argnum + 1;
                pasteup_opts = varargin{argnum};
                args2delete(end+1:end+2) = [argnum-1 argnum];
        end
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];
if strcmp(display_str, 'off') == 1
    varargin{end+1} = 'nodisplay';
end
if trialargnum == 0
    error('lfp_makepasteup:badtrialarg', ...
        'There must be a numeric value or [] supplied for <trials>.' );
end

% Process the pasteup_opts (if any)
argnum = 1;
dataflag = false;
IEIs = [];
if ~isempty(pasteup_opts) && ~iscell(pasteup_opts)
    error('lfp_makepasteup:pasteup_opts', ...
        '''pasteup_opts'' must be followed by a cell array.');
end
while argnum <= length(pasteup_opts)
    switch pasteup_opts{argnum}
        case 'autoselect'
            for k = 1:length(aligns)
                lfp_selectByRule(sprintf( ...
                    'lfp_SelectedTrials(trial) && HasEvent(%s)', ...
                    mat2str(aligns{k}) ));
            end
            args2delete(end+1) = argnum; %#ok<*AGROW>
        case 'aux'
            argnum = argnum + 1;
            auxdata = pasteup_opts{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        case 'blankevtlabels'
            blankevtlabels = true;
        case 'data'
            dataflag = true;
        case 'bg'
            argnum = argnum + 1;
            bg = pasteup_opts{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        case 'intervals'
            argnum = argnum + 1;
            IEIs = pasteup_opts{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        case 'joinmarkers'
            showjoinmarkers = true;
        case 'nodisplay'
            display_str = 'off';
            args2delete(end+1) = argnum;
        case 'title'
            argnum = argnum + 1;
            titlestr = pasteup_opts{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        otherwise
            error('lfp_makepasteup:badoption', ...
                ['The option "' dg_thing2str(pasteup_opts{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isempty(IEIs) && (length(IEIs) ~= length(aligns) - 1)
    error('lfp_makepasteup:IEIs', ...
        'The list of intervals must be one element shorter than the list of alignments');
end

if isempty(trials)
    trials = 1:length(lfp_SelectedTrials);
end
trials = lfp_enabledTrials(trials);
varargin{trialargnum} = trials;

if isempty(lfp_XLimAll) || ~isnumeric(lfp_XLimAll) || lfp_XLimAll(2) < lfp_XLimAll(1) %#ok<*NODEF>
    error('lfp_makepasteup:badxlim', ...
        'Bad value of lfp_XLimAll: %s', dg_thing2str(lfp_XLimAll));
end
old_xlimall = lfp_XLimAll;
old_alignmentref = lfp_AlignmentRef;

try

    if ~dataflag
        hOutFig = figure('Visible', display_str);
        hOutAx = axes;
        hold on;
    end


    if isempty(IEIs)
        % Compute median inter-event intervals
        intervals = zeros(length(trials), length(aligns)-1);
        for trialidx = 1:length(trials)
            trial = trials(trialidx);
            trialevents = lfp_Events( ...
                lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
            thatTS = NaN;
            for alignidx = 1:(length(aligns)-1)
                if isnan(thatTS)
                    thisTS = trialevents(ismember( ...
                        trialevents(:,2), aligns{alignidx} ), 1);
                else
                    if isempty(thatTS)
                        thisTS = [];
                    else
                        thisTS = thatTS(1);
                    end
                end
                thatTS = trialevents( ...
                    ismember(trialevents(:,2), aligns{alignidx+1}), 1 );
                if isempty(thisTS) || isempty(thatTS)
                    if norefOKflag
                        continue
                    else
                        error('lfp_makepasteup:noref3', ...
                            'Trial %d has no alignment event %s.', ...
                            trial, mat2str(aligns{alignidx+1}) );
                    end
                end
                intervals(trialidx, alignidx) = thatTS(1) - thisTS(1);
            end

        end
        % <IEIs2calc> determines the time window of data that goes into the
        % calculation aligned on each of <aligns>.  <offsets> determines the
        % nominal time in seconds at which to plot "time zero" for each of
        % <aligns>, except that the first of <aligns> is by definition plotted at 0
        % seconds.
        IEIs2calc = median(intervals, 1);
        offsets = cumsum(IEIs2calc);
    else
        IEIs2calc = IEIs;
        offsets = cumsum(IEIs2calc);
    end

    % Compute the subplots, pasting results into the pasteup as we go:
    israster = false;
    for alignidx = 1:length(aligns)
        % Set the time interval over which to compute <funch>:
        lfp_AlignmentRef = aligns{alignidx};
        if alignidx > 1
            lfp_XLimAll(1) = -IEIs2calc(alignidx-1)/2;
        end
        if alignidx < length(aligns)
            lfp_XLimAll(2) = IEIs2calc(alignidx)/2;
        else
            lfp_XLimAll(2) = old_xlimall(2);
        end

        if dataflag
            [val, pd] = feval(funch, varargin{:});
            if alignidx == 1
                values = val;
                plotdata = pd;
                plotdata.win = old_xlimall;
                plotdata.offsets = offsets;
                plotdata.aligns = aligns;
                names = fieldnames(values);
                thisoffset = 0;
            else
                thisoffset = offsets(alignidx-1);
                points2use = 1:length(pd.timepts);
                if abs(pd.timepts(1) + thisoffset - plotdata.timepts(end)) ...
                        < lfp_SamplePeriod/2
                    % It's just the same point repeating, so eliminate the
                    % repeated point.
                    points2use(1) = [];
                end
                if pd.timepts(points2use(1)) + thisoffset < ...
                        plotdata.timepts(end)
                    detailstr = sprintf( ...
                        '\nalignidx=%d, align=%s, last=%.6f, next=%.6f', ...
                        alignidx, dg_thing2str(aligns{alignidx}), ...
                        plotdata.timepts(end), ...
                        pd.timepts(1) + thisoffset);
                    warning('lfp_makepasteup:backwards', ...
                        'Time is going backwards!  That''s usually not good...%s', ...
                        detailstr);
                end
                plotdata.timepts = ...
                    [plotdata.timepts pd.timepts(points2use) + thisoffset];
                for nameidx = 1:length(names)
                    values.(names{nameidx}) = ...
                        [ values.(names{nameidx}) ...
                        val.(names{nameidx})(:, points2use) ];
                end
            end
            varargout = {values, plotdata};
        else
            hF(alignidx) = feval(funch, varargin{:});
            % Determine what type of plot it is and check consistency.
            % We identify line plots by process of elimination because any
            % plot that contains event markers will contain line objects
            % for the event markers.
            hData{alignidx} = findobj(hF(alignidx), 'Type', 'image');
            if ~isempty(hData{alignidx})
                thisplottype = 'image';
            else
                hData{alignidx} = findobj(hF(alignidx), 'Type', 'patch');
                if ~isempty(hData{alignidx})
                    thisplottype = 'patch';
                else
                    thisplottype = 'line';
                end
            end
            hA = get(hF(alignidx), 'CurrentAxes');
            if alignidx == 1
                props2copy = ...
                    {'XGrid', 'XMinorGrid', 'YGrid', 'YMinorGrid', 'YDir'};
                for k = 1:length(props2copy)
                    set(hOutAx,  props2copy{k}, get(hA, props2copy{k}));
                end
                texts2copy = {'XLabel', 'YLabel'};
                for k = 1:length(texts2copy)
                    hTnew = copyobj(get(hA, texts2copy{k}), hOutAx);
                    set(hOutAx,  texts2copy{k}, hTnew);
                end
                plottype = thisplottype;
                thisoffset = 0;
            else
                if ~isequal(plottype, thisplottype)
                    error('lfp_makepasteup:plottype', ...
                        'The type of plot returned varies with alignment event' );
                end
                thisoffset = offsets(alignidx-1);
            end

            % Copy relevant data object(s) to pasteup with appropriate x offset.
            hL = findobj(hF(alignidx), 'Type', 'line');
            ylimits = get(hA, 'YLim');
            origxlimits = get(hA, 'XLim');
            % We assume here that every plot may contain an arbitrary number of
            % two-point line objects that represent event markers, and these will
            % appear before any line data in hL.
            switch plottype
                case {'image' 'patch'}
                    if alignidx == 1
                        xlimits = origxlimits;
                    end
                    % patches are delicate little creatures and can't easily be
                    % modified once created; therefore instead of copying and
                    % modifying, we collect all the attributes we need and
                    % create a clonetoid.
                    if isequal(plottype, 'patch')
                        origcolor = get(hData{alignidx}(end), 'FaceColor');
                        if isequal(origcolor, 'flat')
                            % 'flat' does not work as a color
                            origcolor = 'k';
                        end
                        origX = get(hData{alignidx}(end), 'XData');
                        origY = get(hData{alignidx}(end), 'YData');
                        patch(origX + thisoffset, ...
                            origY, origcolor, 'Parent', hOutAx);
                    else
                        hCopy = copyobj(hData{alignidx}(end), hOutAx);
                        set(hCopy, 'XData', ...
                            (get(hCopy, 'XData')) + thisoffset );
                        set(hOutAx, 'YLim', ylimits);
                        if alignidx == length(aligns)
                            xlimits(2) = origxlimits(2) + thisoffset;
                            hCB = findobj(hF(alignidx), 'Tag', 'Colorbar');
                            cbarlabel = get(get(hCB, 'YLabel'), 'String');
                            hCB_new = colorbar('peer', hOutAx);
                            set(get(hCB_new,'YLabel'), 'String', cbarlabel);
                        end
                    end
                    if alignidx > 1
                        if alignidx > 2
                            xposn = (offsets(alignidx-1) + offsets(alignidx-2))/2;
                        else
                            xposn = offsets(alignidx-1)/2;
                        end
                        if showjoinmarkers
                            plot(hOutAx, [xposn xposn], ylimits, 'w', ...
                                'Linewidth', 2);
                        end
                    end
                case 'line'
                    if isempty(hL)
                        error('lfp_makepasteup:nodata', ...
                            'Line plot %d contains no data', alignidx );
                    else
                        % Determine if raster or standard line plot
                        numpts = zeros(length(hL),1);
                        for k = 1:length(hL)
                            numpts(k) = length(get(hL(k), 'XData'));
                        end
                        if ~any(numpts>2)
                            israster = true;
                        end
                        % Draw bg rectangle
                        if israster
                            xdata = origxlimits;
                        else
                            xdata = get(hL(end), 'XData');
                        end
                        if ~isempty(bg)
                            if mod(alignidx, 2)
                                figure(hOutFig);
                                rectangle('Position', ...
                                    [xdata(1) + thisoffset, ylimits(1), ...
                                    xdata(end)-xdata(1), ...
                                    ylimits(2)-ylimits(1) ], ...
                                    'FaceColor', bg{1}, ...
                                    'EdgeColor', 'none' );
                            elseif length(bg) > 1 && mod(alignidx, 2) == 0
                                figure(hOutFig);
                                rectangle('Position', ...
                                    [xdata(1) + thisoffset, ylimits(1), ...
                                    xdata(end)-xdata(1), ...
                                    ylimits(2)-ylimits(1) ], ...
                                    'FaceColor', bg{2}, ...
                                    'EdgeColor', 'none' );
                            end
                        end
                        if ~israster
                            % Copy and delete line data from source, but not event
                            % markers.  (ALL line objects in a raster plot are
                            % treated as event markers.)
                            linedataidx = find(numpts>2);
                            for k = length(linedataidx):-1:1
                                hCopy = copyobj(hL(linedataidx(k)), hOutAx);
                                set(hCopy, 'XData', ...
                                    (get(hCopy, 'XData')) + thisoffset );
                            end
                            hL(linedataidx) = [];
                        end
                    end
                otherwise
                    error('lfp_makepasteup:oops', ...
                        'Oops - internal malfunction' );
            end

            % At this point, hL contains only event markers; copy them if they are
            % within lfp_XLimAll.
            if ~isempty(hL)
                for k = 1:length(hL)
                    xdata = (get(hL(k), 'XData'));
                    if xdata(1) >= lfp_XLimAll(1) && xdata(1) <= lfp_XLimAll(2)
                        hCopy = copyobj(hL(k), hOutAx);
                        set(hCopy, 'XData', ...
                            (get(hCopy, 'XData')) + thisoffset );
                    end
                end
            end
            if strcmp(display_str, 'on') == 1
                figure(hOutFig);
            end
            if isempty(lfp_EventNames{lfp_AlignmentRef(1)}) %#ok<*USENS>
                if blankevtlabels
                    evtname = '';
                else
                    evtname = sprintf('evt%s', dg_thing2str(lfp_AlignmentRef));
                end
            else
                if length(lfp_AlignmentRef) == 1
                    evtname = lfp_EventNames{lfp_AlignmentRef};
                else
                    evtname = sprintf('evt%s', dg_thing2str(lfp_AlignmentRef));
                end
            end
            eventlabels(alignidx) = text(thisoffset, ylimits(2), evtname, ...
                'Parent', hOutAx, ...
                'Interpreter', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom' );
            close(hF(alignidx));
        end
    end

    if ~dataflag
        % Remove the white space from the ends of the plot:
        if isequal(plottype, 'image')
            set(hOutAx, 'XLim', xlimits);
        end

        if ~israster
            % Make sure all event markers extend the full height of the plot;
            % assume all 2-point lines are event markers.
            ylimits = get(hOutAx, 'YLim');
            hL = findobj(hOutAx, 'Type', 'line');
            for k = 1:length(hL)
                ydata = get(hL(k), 'YData');
                if length(ydata) == 2
                    set(hL(k), 'YData', ylimits);
                end
            end
        end

        % Make sure all the bg rectangles extend the full height of the plot
        height = ylimits(2) - ylimits(1);
        hR = findobj(hOutAx, 'Type', 'rectangle');
        for k = 1:length(hR)
            posn = get(hR(k), 'Position');
            posn(2) = ylimits(1);
            posn(4) = height;
            set(hR(k), 'Position', posn);
        end

        % Make sure all the event labels are at the full height of the plot
        for k = 1:length(eventlabels)
            if israster
                yposn = ylimits(1);
            else
                yposn = ylimits(2);
            end
            posn = get(eventlabels(k), 'Position');
            posn(2) = yposn;
            set(eventlabels(k), 'Position', posn);
            set(eventlabels(k), 'Units', 'normalized');
        end

        % Make sure the plot is still the full height of the plot (doh!)  (Note
        % that this disables automatic y scaling)
        set(hOutAx, 'YLim', ylimits);

        set(hOutFig, 'UserData', {varargin auxdata});
        hT = title(hOutAx, sprintf('%s\n', titlestr), 'Interpreter', 'none');
        set(hT, 'ButtonDownFcn', ...
            'lfp_makepasteup_userdata=get(gcf,''UserData''); disp(''UserData in lfp_makepasteup_userdata''); disp(''Arguments:''); disp(dg_thing2str(lfp_makepasteup_userdata{1})); disp(''Aux data:''); disp(dg_thing2str(lfp_makepasteup_userdata{2}));' );
        varargout = {hOutFig, hOutAx};
    end

catch e
    logmsg = sprintf('%s\n%s\n%s', ...
        e.identifier, e.message);
    for stackframe = 1:length(e.stack)
        logmsg = sprintf('%s\n%s\nline %d', ...
            logmsg, e.stack(stackframe).file, ...
            e.stack(stackframe).line);
    end
    disp(logmsg);
end

% Restore old values
lfp_XLimAll = old_xlimall; %#ok<*NASGU>
lfp_AlignmentRef = old_alignmentref;

