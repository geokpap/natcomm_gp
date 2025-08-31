function [hF, hI, hCB, data] = lfp_CSCraster(trials, filenums, window, ...
    varargin)
%[hF, data] = lfp_CSCraster(trials, filenums, window)  Pseudocolor raster
% plot of CSC data
%INPUTS
% trials, filenums, window: as for lfp_disp, etc.  However, only ONE of
%   <trials> of <filenums> may be a non-scalar array.  The one that has
%   multiple elements will be used for the Y axis.  <filenums> may contain
%   elements whose value is NaN, and a row of all NaNs will be displayed at
%   those positions; this is meant to facilitate visually segregating
%   channels into groups.
%OUTPUTS
% hF: figure handle
% hI, hCB: image and colorbar handles as for dg_showGram
% data: the data used to create the image, in trials/filenums x samples
%   format.
%OPTIONS
% 'colors', colors - passes through to dg_compactStripchart when combined
%   with 'stripchart' option.  If <colors> is equal to 'auto', then
%   dg_compactStripchart's default auto-color behavior is invoked.  If
%   'colors' is not given, the default color is middle grey.
% 'filenames' - displays the relevant lfp_FileNames to label each raster
%   when displaying a single trial.
% 'legend', legendstr - when combined with 'stripchart' option, causes a 
%   figure legend to be created, where <legendstr> is a cell string array
%   containing the legend labels for the corresponding line plots.  If
%   <legendstr> is empty, then consecutive integers starting at 1 are used.
%   If <legendstr> is the character string 'names', then the filenames or
%   trial IDs are used.  If 'ovr' is also given, then the lines referenced
%   in the legend comprise the first line of each overlaid set; this
%   implies that if <legendstr> is not empty, it should contain one
%   string for each set of filenums.
% 'nobounds' - ignore the whole issue of trial or recseg boundaries, and
%   just treat all samples as one long continuous recording.  (However: be
%   careful what you wish for.)
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.
% 'offset', offset - <offset> is passed through to dg_compactStripchart
%   when combined with 'stripchart' option.
% 'ovr', morefilenums - <morefilenums> must have as many columns as
%   there are elements in <filenums>.  Each row of <morefilenums> contains
%   one full set of filenums to overlay on the corresponding traces from
%   <filenums> when using 'stripchart' option.  In this case, 'colors' is
%   interpreted as a list of colors one for each full set of filenums in
%   the overlay.
% 'stripchart' - instead of displaying a pseudocolor raster image, displays
%   a dg_compactStripchart.  Each trace in the stripchart can be clicked to
%   display its filename or trial ID in the command window.
%NOTES
% The Y axis is labeled "Trials" or "Channels" as appropriate, but the Y
%   values displayed by the data picker are actually the indices into
%   <trials> or <filenums>, so they will only be the same if <trials> or
%   <filenums> is a consecutive sequence starting at 1.

%$Rev: 361 $
%$Date: 2015-08-10 15:58:30 -0400 (Mon, 10 Aug 2015) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin <3
    window = [];
elseif ~(isa(window, 'double') ...
        && isequal(size(window), [1 2]) ) ...
        && ~isempty(window)
    error('lfp_CSCraster:badwindow', ...
        '<window> must be 1x2 number array.' );
end
if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
if ischar(trials)
    trials = lfp_parseTrialStr(trials);
end
% <trials> is now numeric, i.e. trialnums.
if any(trials > length(lfp_SelectedTrials))
    warning('lfp_CSCraster:trials', ...
        'Ignoring trials %s, which are beyond the last trial.', ...
        dg_canonicalSeries(trials(trials > length(lfp_SelectedTrials))) );
    trials(trials > length(lfp_SelectedTrials)) = [];
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

dg_compactStripchart_opts = {};
defaultcolor = [.5 .5 .5];
colors = defaultcolor;
dBflag = false;
display_str = 'on';
filenamesflag = false;
legendflag = false;
legendstr = {};
lfp_findCommonTime_opts = {'recseg'};
markersize = 12;
morefilenums = [];
normflag = false;
offset = [];
stripchartflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'colors'
            argnum = argnum + 1;
            if argnum > length(varargin) ...
                    || ~isempty(varargin{argnum}) ...
                    && size(varargin{argnum}, 2) ~= 3 ...
                    && ~isequal(varargin{argnum}, 'auto')
                error('lfp_CSCraster:badcolors', ...
                    '''colors'' option requires a 3-column numeric array or ''auto''');
            end
            colors = varargin{argnum};
        case 'dB'
            dBflag = true;
        case 'filenames'
            filenamesflag = true;
        case 'nobounds'
            lfp_findCommonTime_opts = {'nobounds'};
        case 'nodisplay'
            display_str = 'off';
        case 'legend'
            argnum = argnum + 1;
            legendflag = true;
            if argnum > length(varargin) || ...
                    ~isempty(varargin{argnum}) ...
                    && ~iscell(varargin{argnum}) ...
                    && ~isequal(varargin{argnum}, 'names')
                error('lfp_CSCraster:badlegend', ...
                    '''legend'' option requires an argument (see header comments)');
            end
            legendstr = varargin{argnum};
        case 'norm'
            normflag = true;
        case 'offset'
            argnum = argnum + 1;
            offset = varargin{argnum};
        case 'ovr'
            argnum = argnum + 1;
            morefilenums = varargin{argnum};
            hL2 = NaN(size(morefilenums,1) + 1, 1);
        case 'stripchart'
            stripchartflag = true;
        otherwise
            error('lfp_CSCraster:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isequal(colors, 'auto') 
    if isempty(morefilenums)
        dg_compactStripchart_opts(1, end+1:end+2) = {'colors', colors};
    else
        if size(colors,1) > 1
            dg_compactStripchart_opts(1, end+1:end+2) = ...
                {'colors', colors(1,:)};
        else
            dg_compactStripchart_opts(1, end+1:end+2) = {'colors', colors};
        end
    end
end

if numel(trials) > 1 && numel(filenums) > 1
    error('lfp_CSCraster:noscalar', ...
        'At least one of <trials> or <filenums> must be a scalar.');
end

if ~isempty(morefilenums) && size(morefilenums,2) ~= length(filenums)
    error('lfp_CSCraster:morefilenums', ...
        'Number of columns in <morefilenums> must match<filenums>.');
end

[data, interval, trialevents, reftime, window] = lfp_CSCraster_gatherdata( ...
    filenums, trials, window, lfp_findCommonTime_opts);
hF = figure('Visible', display_str);
hA = axes('Parent', hF);
if numel(trials) > 1
    numrows = length(trials);
else
    numrows = length(filenums);
end
if normflag
    data = data / max(data(:));
end
if dBflag
    data = 10 * log10(data);
end
if stripchartflag
    [hF, hA, hL, offset] = dg_compactStripchart(data, offset, ...
        'axes', hA, 'samplesize', lfp_SamplePeriod, ...
        'x0', -interval(1) + 1, dg_compactStripchart_opts{:}); 
    if ~isempty(morefilenums)
        hL2(1) = hL(1);
    end
    rowheight = -offset;
    firstrowheight = offset;
    for k = 1:length(hL)
        if numel(trials) > 1
            detailstr = lfp_getTrialID(trials(k));
        else
            detailstr = lfp_FileNames{filenums(k)}; %#ok<*USENS>
        end
        set(hL(k), ...
            'ButtonDownFcn', ...
            ['fprintf(1,''' detailstr '\n'')'] );
    end
    set(hA, 'YLim', ...
        [min(data(end,:))-(size(data,1)-1)*offset max(data(1,:))]);
    if numel(trials) > 1
        set(hA, 'YTick', []);
    else
        set( hA, 'YTick', ...
            (length(filenums):-1:1) * rowheight + firstrowheight );
        if filenamesflag
            yticklabels = cell(size(filenums));
            for k=1:length(filenums)
                if isnan(filenums(length(filenums) - k + 1))
                    yticklabels{k} = '';
                else
                    yticklabels{k} = lfp_FileNames{ filenums( ...
                        length(filenums) - k + 1) }; 
                end
            end
        else
            yticklabels = filenums(end:-1:1);
        end
        set(hA, 'YTickLabel', yticklabels);
    end
    if ~isempty(morefilenums)
        for rownum = 1:size(morefilenums,1)
            [moredata, interval, trialevents, reftime, window] = ...
                lfp_CSCraster_gatherdata( ...
                morefilenums(rownum,:), trials, window, ...
                lfp_findCommonTime_opts); 
            if size(colors,1) > 1
                if rownum + 1 > size(colors,1)
                    more_Stripchart_opts = {'colors', defaultcolor};
                else
                    more_Stripchart_opts = {'colors', colors(rownum+1,:)};
                end
            else
                more_Stripchart_opts = {'colors', colors * 0.5^rownum};
            end
            [hF3, hA3, hL3] = dg_compactStripchart(moredata, offset, ...
                'axes', hA, 'samplesize', lfp_SamplePeriod, ...
                'x0', -interval(1) + 1, more_Stripchart_opts{:}); %#ok<ASGLU>
            hL2(rownum+1) = hL3(1);
            for k = 1:length(hL3)
                if numel(trials) > 1
                    detailstr = lfp_getTrialID(trials(k));
                else
                    detailstr = lfp_FileNames{filenums(k)}; %#ok<*USENS>
                end
                set(hL3(k), ...
                    'ButtonDownFcn', ...
                    ['fprintf(1,''' detailstr '\n'')'] );
            end
        end
    end
    if legendflag
        if isempty(morefilenums)
            % Normal (non-overlaid) plot legend
            if isempty(legendstr)
                for k = 1:numrows
                    legendstr{k} = int2str(k); %#ok<*AGROW>
                end
            elseif isequal(legendstr, 'names')
                legendstr = {};
                if numel(trials) > 1
                    % Use trial IDs
                    for k = 1:numrows
                        legendstr{k} = lfp_getTrialID(trials(k));
                    end
                else
                    % Use file names
                    for k = 1:numrows
                        legendstr{k} = lfp_FileNames{filenums(k)};
                    end
                end
            else
                if length(legendstr) > numrows
                    legendstr(numrows+1:end) = [];
                else
                    legendstr(end+1:numrows) = {''};
                end
            end
            hleg = legend(hA, legendstr);
        else
            % Overlaid plot; only use first channel from each set.
            numsets = size(morefilenums,1) + 1;
            if isempty(legendstr)
                for k = 1 : numsets
                    legendstr{k} = int2str(k); %#ok<*AGROW>
                end
            elseif isequal(legendstr, 'names')
                legendstr = lfp_FileNames(filenums(1));
                % Use file names
                for k = 2: numsets
                    legendstr{k} = lfp_FileNames{morefilenums(k-1,1)};
                end
            else
                if length(legendstr) > numsets
                    legendstr(numsets+1:end) = [];
                else
                    legendstr(end+1:numsets) = {''};
                end
            end
            hleg = legend(hA, hL2, legendstr);
        end
        set(hleg,'Interpreter','none');
    end
    [mantissa, exp] = dg_findScale(offset);
    calibsize = mantissa * 10^exp;
    calibstr = sprintf(' %1.0e %s', calibsize, ...
        lfp_SamplesUnits{filenums(1)});
    myxlim = get(hA, 'XLim');
    calibX = (myxlim(2)-myxlim(1))/100 + myxlim(1);
    calibY = (calibsize - offset)/2 + [-calibsize 0];
    plot(hA, [calibX calibX], calibY, 'k', 'LineWidth', 1.5);
    text(calibX, mean(calibY), calibstr, 'Parent', hA);
else
    rowheight = 1;
    firstrowheight = 0;
    tvals = lfp_SamplePeriod * ((0:size(data,2)-1) + interval(1));
    hI = imagesc(tvals, [1 numrows], data);
    set(hA, 'YDir', 'reverse', 'YTick', []);
    hCB = colorbar;
    if ~isempty(lfp_CLimAll)
        hCB = dg_recolorGram(hCB, lfp_CLimAll, hI);
    end
    if numel(trials) > 1
        set(get(hCB,'YLabel'), 'String', sprintf( ...
            '%s, %s', lfp_FileNames{filenums}, lfp_SamplesUnits{filenums} )); 
    else
        CBlabel = 'CSC values';
        if normflag
            CBlabel = [CBlabel ' normalized'];
        end
        if dBflag
            CBlabel = [CBlabel ', dB'];
        end
        set(get(hCB,'YLabel'), 'String', CBlabel);
        set(hA, 'YTick', 1:length(filenums));
        if filenamesflag
            yticklabels = cell(size(filenums));
            for k=1:length(filenums)
                if isnan(filenums(k))
                    yticklabels{k} = '';
                else
                    yticklabels{k} = lfp_FileNames{filenums(k)};
                end
            end
        else
            yticklabels = filenums;
        end
        set(hA, 'YTickLabel', yticklabels);
    end
    set(get(hCB,'YLabel'), 'Interpreter', 'none');
end
if numel(trials) > 1
    set(get(hA,'YLabel'), 'String', sprintf( ...
        '%s\nTrials', lfp_FileNames{filenums}));
else
    set(get(hA,'YLabel'), 'String', 'Channels');
end
set(get(hA,'XLabel'), 'String', 'Time, seconds');
set(get(hA,'YLabel'), 'Interpreter', 'none');
set(get(hA,'YLabel'), 'Interpreter', 'none');
hold on;

if numel(trials) == 1
    lfp_plotEvtMarkers( hA, ...
        trialevents, ...
        'reftime', reftime );
else
    % Plot trial-by-trial event markers
    for trialidx = 1:length(trials)
        trial = trials(trialidx);
        trialevents = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), :); %#ok<NODEF>
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), ...
            1 );
        if length(reftime) > 1
            warning('lfp_CSCraster:multiref', ...
                'Trial %d has more than one reference event', trial);
            reftime = reftime(1);
        end
        if isempty(window)
            starttime = lfp_Events( ...
                lfp_TrialIndex(trial,1), 1 );
            endtime = lfp_Events( ...
                lfp_TrialIndex(trial,2), 1 );
        else
            starttime = max( ...
                reftime + window(1), ...
                lfp_Events( ...
                lfp_TrialIndex(trial,1), 1 ));
            endtime = min( ...
                reftime + window(2), ...
                lfp_Events( ...
                lfp_TrialIndex(trial,2), 1 ));
        end
        eventrange = reshape(find( trialevents(:,1) >= starttime ...
            & trialevents(:,1) <= endtime ), 1, []);
        for myevtix = eventrange
            evtid = trialevents(myevtix,2);
            if lfp_SelectedEventIDs(evtid)
                eventname = ''; % required for detailstr
                if evtid <= length(lfp_EventNames) 
                    eventname = lfp_EventNames{evtid};
                end
                if evtid <= length(lfp_EventColors) ...
                        && ~isempty(lfp_EventColors{evtid}) 
                    eventcolor = lfp_EventColors{evtid};
                else
                    eventcolor = lfp_EventDefaultColor;
                end
                if evtid <= length(lfp_EventShapes) ...
                        && ~isempty(lfp_EventShapes{evtid}) 
                    eventshape = lfp_EventShapes{evtid};
                else
                    eventshape = lfp_EventDefaultShape;
                end
                hL = plot(hA, ...
                    trialevents(myevtix,1) - reftime, ...
                    trialidx * rowheight + firstrowheight, ...
                    'Color', eventcolor, ...
                    'Marker', eventshape, 'MarkerSize', markersize );
                if isequal(class(eventcolor), 'char')
                    eventcolorstr = eventcolor;
                else
                    eventcolorstr = mat2str(eventcolor);
                end
                detailstr = sprintf( ...
                    '\\nTimestamp=%.6f\\nLineColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
                    trialevents(myevtix,1), ...
                    eventcolorstr, eventname, evtid, evtid );
                set(hL, ...
                    'ButtonDownFcn', ...
                    ['fprintf(1,''' detailstr '\n'')'] );
            end
        end
    end
end

lfp_createFigTitle(hA, 'CSCras', trials, window, '', ''); 

end


function [data, interval, trialevents, reftime, window] = ...
    lfp_CSCraster_gatherdata(filenums, trials, window, ...
    lfp_findCommonTime_opts)
global lfp_Samples lfp_XLimAll lfp_TrialIndex lfp_Events lfp_AlignmentRef
global lfp_SamplePeriod
trialevents = [];
reftime = [];

if isempty(lfp_XLimAll) && isempty(window)
    % Use whole trials and pad short ones with NaNs as needed
    interval = [];
    rawtrialinfo = [];
    for trial = reshape(trials, 1, [])
        trialevents = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), ...
            : ); 
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), 1 );
        if isempty(reftime)
            refpoint = 0;
        else
            refpoint = lfp_time2index(reftime(1));
        end
        rawtrialinfo(end+1,:) = ...
            [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint]; 
    end
else
    if isempty(window)
        window = lfp_XLimAll;
    end
    [interval, rawtrialinfo] = lfp_findCommonTime(trials, ...
        lfp_findCommonTime_opts{:});
    xlimpoints = round(window/lfp_SamplePeriod);
    interval(1) = max(xlimpoints(1), interval(1));
    interval(2) = min(xlimpoints(2), interval(2));
    if numel(trials) == 1
        trialevents = lfp_Events( ...
            lfp_TrialIndex(trials,1) : lfp_TrialIndex(trials,2), ...
            : ); 
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), 1 );
        trialevents = lfp_Events( ...
            lfp_Events(:,1) >= reftime + window(1) ...
            & lfp_Events(:,1) <= reftime + window(2), ...
            : );
    end
end
if any(rawtrialinfo(:,3)==0)
    error('lfp_CSCraster:noref', ...
        'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(rawtrialinfo(:,3)==0)) );
end
if numel(trials) > 1
    numrows = length(trials);
else
    numrows = length(filenums);
end
if isempty(interval) || numel(filenums) > 1
    % For whole trials or multiple channels, construct the <data> array row
    % by row.
    if isempty(interval)
        ptsBeforeRef = rawtrialinfo(:,3) - rawtrialinfo(:,1);
        ptsAfterRef = rawtrialinfo(:,2) - rawtrialinfo(:,3);
    else
        ptsBeforeRef = min( ...
            rawtrialinfo(:,3) - rawtrialinfo(:,1), -interval(1) );
        ptsAfterRef = min( ...
            rawtrialinfo(:,2) - rawtrialinfo(:,3), interval(2) );
    end        
    numpts = max(ptsBeforeRef) + max(ptsAfterRef) + 1;
    data = NaN(numrows, numpts);
    refpt = max(ptsBeforeRef) + 1;
    startpt = refpt - ptsBeforeRef;
    endpt = refpt + ptsAfterRef;
    if numel(trials) > 1
        for trialidx = 1:length(trials)
            data(trialidx, startpt(trialidx):endpt(trialidx)) = ...
                lfp_Samples{filenums}(rawtrialinfo(trialidx,1) : ...
                rawtrialinfo(trialidx,2)); 
        end
    else
        for filenumidx = 1:length(filenums)
            if isnan(filenums(filenumidx))
                data(filenumidx, startpt:endpt) = NaN;
            else
                data(filenumidx, startpt:endpt) = lfp_Samples{ ...
                    filenums(filenumidx) }( ...
                    rawtrialinfo(3) - ptsBeforeRef : ...
                    rawtrialinfo(3) + ptsAfterRef );
            end
        end
    end
    interval = [-max(ptsBeforeRef) max(ptsAfterRef)];
else
    % This should run faster for large numbers of samples
    idxrange = interval(1):interval(2);
    indices = (repmat(idxrange, size(rawtrialinfo(:,3))) ...
        + repmat(rawtrialinfo(:,3), size(idxrange)) );
    if isempty(indices)
        error('lfp_CSCraster:nodata', ...
            ['No samples were selected; note that if lfp_XLimAll is\n' ...
            'empty, <window> is clipped to start and end of trial.'])
    end
    data = reshape(lfp_Samples{filenums}(indices), size(indices)); 
end
end
