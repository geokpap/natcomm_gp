function hF = lfp_disp2d(varargin)
%hF = lfp_disp2d(trials, filenums, window)
% Makes plots of 2-D tracker data (e.g. eye tracker or Neuralynx video
% tracker).
%INPUTS
% trials, window: as usual.
% filenums: a two-element or 3-element vector of filenums to use as data
%   sources, filenums(1) for X, filenums(2) for Y, filenums(3) optionally
%   for color coding of each point plotted (NaN values are not plotted in
%   this case.)
%OUTPUTS
% hF: handle to newly plotted figure.
%OPTIONS
% All options that can be given to lfp_getSamples can be given to
% lfp_disp2d with the same effect.  In addition, the following options are
% available:
% 'axes', hA - plots into the pre-existing axes handle <hA>.
% 'bg', bgargs - before plotting data, draws a set of background graphics
%   objects as specified by <bgargs>, which is a cell vector. Each cell is
%   itself a cell vector that starts with the name of a Matlab function to
%   invoke to draw something in the same axes where the data will be
%   plotted.  The function specified must accept the property-value pair
%   "'Parent', hA" as part of its argument list (e.g. the 'patch' or
%   'rectangle' functions).  The remainder of the cell vector contains
%   arguments to the function, not including "'Parent', hA" (which is
%   automatically added at the end). Example:
%       'bg', {{'rectangle'} {'patch', [.1 .2 .3], [.1 .3 .1], 'r'}}
%   draws a default rectangle and then a red-filled triangle with corners
%   at (.1, .1), (.2, .3), (.3, .1).
% 'colorgrid', gridspacing - when using 'colorTS' or a third filenum to
%   color-code the points, this option divides the X-Y plane into squares
%   of width <gridspacing> and colors each square according to the average
%   value for all points whose (X, Y) coordinates fall within the square.
%   Boundaries are handled as in Matlab's 'histc' function, i.e. using the
%   comparison operators >= and <.  If a bin contains more than half NaN
%   values, a warning is given.
% 'events' - plot event markers on top of the track at the positions
%   interpolated between the bracketing samples.
% 'eventsmarker', marker - plot events using the specified marker symbol
%   (see Matlab Line Properties); default is '^'.
% 'colorTS' - code the timestamp of each point as a color.
% 'marker', marker - plot using the specified marker symbol (see Matlab
%   'scatter' function or Line Properties); default is 'o'.
% 'markersize', S - use markers of size <S> (see Matlab 'scatter' function;
%   when plotting without color coding, 'plot' is used instead and sqrt(S)
%   is used for 'MarkerSize').
% 'minpts', N - only applies when used together with 'colorgrid'; sets the
%   minimum number of points in a grid cell that is required to include it
%   in the plot.  Default is 3.  When 'trialweight' is used, each trial
%   counts as one point.
% 'movie', avifilename - create a Matlab movie using dg_plotmovie and save
%   it as a 25 frames/sec (same rate as PAL/SECAM) AVI file at the location
%   specified by <avifilename>, which can be an absolute or relative path
%   to the output file.  Heedlessly blows away any pre-existing file at
%   that location.  Does not work with 'colorgrid' or "whole-trial-data
%   mode" (i.e. 'notrunc' with empty window and lfp_XLimAll).  Do not touch
%   the figure window that is creating the movie until it is finished (see
%   dg_plotmovie).
% 'movtimescale', scale - magnifies the time scale of a movie made using
%   the 'movie' option by a factor of <scale> (e.g. scale = 2 means that
%   the movie will take twice as long to play as it took to record the
%   data).  <scale> can be any positive floating-point number.
% 'subsample', N - plot every Nth point instead of all points.
% 'trialweight' - only applies when used together with 'colorgrid'.  
%   Computes the average in each bin by weighting each trial equally,
%   rather than the default computation which weights each point equally.
%   The distinction is that if the track stays in one bin for an unusually
%   long time on one trial, then by default that trial contributes
%   proportionally more to the average than the other trials; but when
%   using 'trialweight', the average value of all the points in the bin for
%   a given trial contributes equally to the final average, regardless of
%   how many points there were in the bin on that trial.
%NOTES
% If all you need is to get the data plotted in <hF>, the best way is to
% call lfp_getSamples directly.

%$Rev: 396 $
%$Date: 2019-02-15 18:13:19 -0500 (Fri, 15 Feb 2019) $
%$Author: dgibson $

global lfp_FileNames lfp_SamplePeriod lfp_XLimAll lfp_Samples ...
    lfp_AlignmentRef lfp_Events lfp_TrialIndex

[trials, filenums, window, arglist, verboseflag, ...
    getSamplesOpts] = lfp_CSCboilerplate(varargin);
[sampledata, timepts, ~, evtidx, badtrials, trials, filenums, triginfo] = ...
    lfp_getSamples(trials, filenums, window, getSamplesOpts{:});
trials = setdiff(trials, badtrials);
if isempty(getSamplesOpts) && isempty(arglist)
    optstring = '';
else
    optstring = dg_thing2str([getSamplesOpts arglist]);
end

argnum = 1;
avifilename = '';
bgargs = {};
colorTSflag = false;
gridspacing = 0;
eventsflag = false;
eventsmarker = '^';
marker = 'o';
markersize = 12;
trialweightflag = false;
N = 3;
hA = [];
Nsub = 1;
scale = 1;
while argnum <= length(arglist)
    if ischar(arglist{argnum})
        switch arglist{argnum}
            case 'axes'
                argnum = argnum + 1;
                hA = arglist{argnum};
            case 'bg'
                argnum = argnum + 1;
                bgargs = arglist{argnum};
            case 'colorgrid'
                argnum = argnum + 1;
                gridspacing = arglist{argnum};
            case 'events'
                eventsflag = true;
            case 'eventsmarker'
                argnum = argnum + 1;
                eventsmarker = arglist{argnum};
            case 'colorTS'
                colorTSflag = true;
            case 'marker'
                argnum = argnum + 1;
                marker = arglist{argnum};
            case 'markersize'
                argnum = argnum + 1;
                markersize = arglist{argnum};
            case 'minpts'
                argnum = argnum + 1;
                N = arglist{argnum};
                if ~isnumeric(N) || ~isscalar(N)
                    error('lfp_disp2d:N', ...
                        'The value following ''minpts'' must be a numeric scalar.');
                end
            case 'movie'
                argnum = argnum + 1;
                avifilename = arglist{argnum};
            case 'movtimescale'
                argnum = argnum + 1;
                scale = arglist{argnum};
            case 'subsample'
                argnum = argnum + 1;
                Nsub = arglist{argnum};
            case 'trialweight'
                trialweightflag = true;
            otherwise
                error('lfp_disp2d:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(arglist{argnum}));
        end
    else
        error('lfp_disp2d:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(arglist{argnum}));
    end
    argnum = argnum + 1;
end

if numel(filenums) < 2 || numel(filenums) > 3
    error('lfp_disp2d:filenums', ...
        '<filenums> must be a 2- or 3-element vector.');
end
if colorTSflag && length(filenums) == 3
    error('lfp_disp2d:color', ...
        'The options ''colorCSC'' and ''colorTS'' are mutually exclusive.');
end
if gridspacing ~= 0 
    if ~colorTSflag && length(filenums) ~= 3
        error('lfp_disp2d:gridspacing', ...
            'The ''gridspacing'' option has no effect unless ''colorCSC'' or ''colorTS'' is used.');
    elseif gridspacing ~= 0 && ~isempty(avifilename)
        error('lfp_disp2d:gridspacing', ...
            'The ''gridspacing'' and ''movie'' options are mutually exclusive.');
    end
end
if ~isempty(avifilename)
    if gridspacing ~= 0
        error('lfp_disp2d:optsconf1', ...
            '''movie'' does not work with ''colorgrid''.');
    end
    if any(cellfun(@isequal, getSamplesOpts, ...
            repmat({'notrunc'}, size(getSamplesOpts)) )) && ...
            isempty(window) && isempty(lfp_XLimAll)
        error('lfp_disp2d:optsconf2', ...
            ['''movie'' does not work with "whole-trial-data mode" ' ...
            '(i.e.\nwith ''notrunc'' when <window> and <lfp_XLimAll> ' ...
            'are both empty).']);
    end
end
if eventsflag && length(trials) ~= size(triginfo, 1)
    warning('lfp_disp2d:events', ...
            '''events'' does not work with ''multitrig''.');
        eventsflag = false;
end

% find xlim and ylim
if iscell(sampledata)
    xmin = Inf;
    xmax = -Inf;
    ymin = Inf;
    ymax = -Inf;
    for tridx = 1:length(sampledata)
        xmin = min([xmin
            reshape(sampledata{tridx}(1:Nsub:end, 1), [], 1) ]);
        ymin = min([ymin
            reshape(sampledata{tridx}(1:Nsub:end, 2), [], 1) ]);
        xmax = max([xmax
            reshape(sampledata{tridx}(1:Nsub:end, 1), [], 1) ]);
        ymax = max([ymax
            reshape(sampledata{tridx}(1:Nsub:end, 2), [], 1) ]);
    end
else
    for trignum = 1:size(sampledata,2)
        xmin = min(min(sampledata(1:Nsub:end, :, 1)));
        xmax = max(max(sampledata(1:Nsub:end, :, 1)));
        ymin = min(min(sampledata(1:Nsub:end, :, 2)));
        ymax = max(max(sampledata(1:Nsub:end, :, 2)));
    end
end

if gridspacing ~= 0
    % Set up the grid.  The "vals" are to be used as bin edges.
    xvals = ( xmin/gridspacing + ...
        (0 : (1 + (xmax - xmin) / gridspacing)) ) * gridspacing;
    yvals = ( ymin/gridspacing + ...
        (0 : (1 + (ymax - ymin) / gridspacing)) ) * gridspacing;
    numx = length(xvals) - 1;
    numy = length(yvals) - 1;
    colorgrid = cell(numx, numy);
end

if isempty(hA)
    hF = figure;
    hA = axes('Parent', hF);
    axis(hA, 'equal');
else
    hF = get(hA, 'Parent');
end
set(hA, 'NextPlot', 'add', 'YDir', 'reverse');
xlabel(hA, lfp_FileNames{filenums(1)}, 'Interpreter', 'none');
ylabel(hA, lfp_FileNames{filenums(2)}, 'Interpreter', 'none');
lfp_createFigTitle(hA, '2D Tracker', trials, window, optstring, '');
if ~isempty(bgargs)
    for obnum = 1:length(bgargs)
        feval(bgargs{obnum}{1}, bgargs{obnum}{2:end}, 'Parent', hA);
    end
end
if isempty(get(hA, 'Children'))
    set(hA, 'XLim', [xmin xmax]);
    set(hA, 'YLim', [ymin ymax]);
end

mov = [];
fps = 25;
if ~isempty(avifilename)
    framesize = round(1 / (Nsub * scale * lfp_SamplePeriod * fps));
end
if verboseflag
    fprintf('Creating plot...\n');
end
if colorTSflag || length(filenums) == 3
    % Set up color display parameters, graphics handles, etc.
    hCB = colorbar('peer', hA);
    cmin = Inf;
    cmax = -Inf;
    if length(filenums) == 3
        % third channel codes color
        if iscell(sampledata)
            % whole-trial-data mode
            for tridx = 1:length(timepts)
                cmin = min([ min(min(sampledata{tridx}(1:Nsub:end, 3))) ...
                    cmin ]);
                cmax = max([ max(max(sampledata{tridx}(1:Nsub:end, 3)))...
                    cmax ]);
            end
            caxis(hA, [cmin cmax]);
        else
            caxis(hA, [ min(min(sampledata(1:Nsub:end, :, 3))) ...
                max(max(sampledata(1:Nsub:end, :, 3))) ]);
        end
        ylabel(hCB, lfp_FileNames{filenums(3)}, 'Interpreter', 'none');
    else
        % color timestamps
        if iscell(sampledata)
            % whole-trial-data mode
            for tridx = 1:length(timepts)
                cmin = min([ min(min(timepts{tridx}(1:Nsub:end))) ...
                    cmin ]);
                cmax = max([ max(max(timepts{tridx}(1:Nsub:end)))...
                    cmax ]);
            end
        else
            caxis(hA, [ min(timepts) max(timepts) ]);
        end
        ylabel(hCB, 'Time, s', 'Interpreter', 'none');
    end
    % Gather color data for each trigger
    if colorTSflag && ~iscell(sampledata)
        % There is only one time scale that is common to all trials:
        colordata = timepts(1:Nsub:end);
    end
    for trignum = 1:size(sampledata,2)
        if colorTSflag && iscell(sampledata)
            % Each trial has its own timescale:
            colordata = timepts{trignum}(1:Nsub:end); 
        end
        if length(filenums) == 3
            if iscell(sampledata)
                % whole-trial-data mode
                colordata = sampledata{trignum}( ...
                    1:Nsub:end, 3 );
            else
                colordata = sampledata(1:Nsub:end, trignum, 3);
            end
        end
        if gridspacing == 0
            % Normal default plot
            if isempty(avifilename)
                if iscell(sampledata)
                    % whole-trial-data mode
                    scatter(hA, ...
                        sampledata{trignum}(1:Nsub:end, 1), ...
                        sampledata{trignum}(1:Nsub:end, 2), ...
                        markersize, colordata, 'fill', marker);
                else
                    scatter(hA, sampledata(1:Nsub:end, trignum, 1), ...
                        sampledata(1:Nsub:end, trignum, 2), ...
                        markersize, colordata, 'fill', marker);
                end
            else
                % 'movie' has been invoked
                if isempty(mov)
                    % Create new movie on first trigger:
                    mov = dg_plotmovie(framesize, ...
                        sampledata(1:Nsub:end, trignum, 1), ...
                        sampledata(1:Nsub:end, trignum, 2), ...
                        'scatter', markersize, colordata, ...
                        'fill', marker, 'axes', hA);
                else
                    % Append to old movie on subsequent triggers:
                    mov = dg_plotmovie(framesize, ...
                        sampledata(1:Nsub:end, trignum, 1), ...
                        sampledata(1:Nsub:end, trignum, 2), ...
                        'scatter', markersize, colordata, ...
                        'fill', marker, 'axes', hA, ...
                        'append', mov);
                end
            end
        else
            % 'colorgrid' has been invoked. Instead of plotting at this
            % point, we just gather data across triggers here:
            for xidx = 1:numx
                for yidx = 1:numy
                    if iscell(sampledata)
                        % whole-trial-data mode
                        isinbin = ...
                            sampledata{trignum}(1:Nsub:end, 1) >= xvals(xidx) ...
                            & sampledata{trignum}(1:Nsub:end, 1) < xvals(xidx+1) ...
                            & sampledata{trignum}(1:Nsub:end, 2) >= yvals(yidx) ...
                            & sampledata{trignum}(1:Nsub:end, 2) < yvals(yidx+1);
                    else
                        isinbin = ...
                            sampledata(1:Nsub:end, trignum, 1) >= xvals(xidx) ...
                            & sampledata(1:Nsub:end, trignum, 1) < xvals(xidx+1) ...
                            & sampledata(1:Nsub:end, trignum, 2) >= yvals(yidx) ...
                            & sampledata(1:Nsub:end, trignum, 2) < yvals(yidx+1);
                    end
                    if trialweightflag
                        binmean = nanmean(colordata(isinbin));
                        if ~isnan(binmean)
                            colorgrid{xidx, yidx} = [ colorgrid{xidx, yidx}
                                binmean ];
                        end
                    else
                        colorgrid{xidx, yidx} = [ colorgrid{xidx, yidx}
                            reshape(colordata(isinbin), [], 1) ];
                    end
                end
            end
        end
        if eventsflag
            trialevtidx = evtidx{trignum};
            if iscell(sampledata)
                % whole-trial-data mode, plot all events
            else
                % remove events outside of analysis window
                trialnum = trials(trignum);
                relalignidx = find(ismember(lfp_Events( ...
                    lfp_TrialIndex(trialnum,1):lfp_TrialIndex(trialnum,2), ...
                    2 ), lfp_AlignmentRef));
                alignTS = lfp_Events( ...
                    relalignidx + lfp_TrialIndex(trialnum,1) - 1, 1 );
                evts2del = ...
                    lfp_Events(trialevtidx, 1) < alignTS + window(1) ...
                    | lfp_Events(trialevtidx, 1) >= alignTS + window(2);
                trialevtidx(evts2del) = [];
            end
            % plot event markers at location interpolated between samples
            for trialevtidx2 = 1:length(trialevtidx)
                eventTS = lfp_Events(trialevtidx(trialevtidx2), 1);
                eventID = lfp_Events(trialevtidx(trialevtidx2), 2);
                eventsampidx = lfp_time2index(eventTS);
                if lfp_index2time(eventsampidx) <= eventTS
                    % eventsampidx comes before or at event time
                    samp1 = eventsampidx;
                    samp2 = eventsampidx + 1;
                else
                    samp1 = eventsampidx - 1;
                    samp2 = eventsampidx;
                end
                TS1 = lfp_index2time(samp1);
                TS2 = lfp_index2time(samp2);
                [evtname, evtcolor, evtcolorstr] = ...
                    lfp_getEvtProps(eventID, {});
                detailstr = sprintf( ...
                    '\\nTimestamp=%.6f\\nMarkerColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
                    eventTS, evtcolorstr, evtname, eventID, eventID );
                plot( interp1( [TS1 TS2], ...
                    [ lfp_Samples{filenums(1)}(samp1) ...
                    lfp_Samples{filenums(1)}(samp2) ], ...
                    eventTS ), ...
                    interp1( [TS1 TS2], ...
                    [ lfp_Samples{filenums(2)}(samp1) ...
                    lfp_Samples{filenums(2)}(samp2) ], ...
                    eventTS ), ...
                    'ButtonDownFcn', ...
                    ['fprintf(1,''' detailstr '\n'')'], ...
                    'Marker', eventsmarker, ...
                    'MarkerFaceColor',  evtcolor, ...
                    'MarkerEdgeColor',  evtcolor);
            end
        end
    end
else
    % No color coding
    if isempty(avifilename)
        if iscell(sampledata)
            % whole-trial-data mode
            for tridx = 1:length(sampledata)
                plot(hA, reshape(sampledata{tridx}(1:Nsub:end,1), [], 1), ...
                    reshape(sampledata{tridx}(1:Nsub:end,2), [], 1), ...
                    'LineStyle', 'none', 'Marker', marker, ...
                    'MarkerSize', sqrt(markersize), ...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            end
        else
            plot(hA, reshape(sampledata(1:Nsub:end,:,1), [], 1), ...
                reshape(sampledata(1:Nsub:end,:,2), [], 1), ...
                'LineStyle', 'none', 'Marker', marker, ...
                'MarkerSize', sqrt(markersize), ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        end
    else
        % 'movie' has been invoked
        mov = dg_plotmovie(framesize, ...
            reshape(sampledata(1:Nsub:end,:,1), [], 1), ...
            reshape(sampledata(1:Nsub:end,:,2), [], 1), ...
            'LineStyle', 'none', 'Marker', marker, ...
            'MarkerSize', sqrt(markersize), ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'axes', hA);
    end
end
if gridspacing == 0
    % Convert from Matlab "movie" to AVI file:
    if verboseflag
        fprintf('Creating movie file %s\n', avifilename);
    end
    if ~isempty(mov)
        v = ver('matlab');
        vtok = regexp(v.Version, '^(\d+\.\d+)', 'tokens');
        if str2double(vtok{1}) >= 7.12
            vObj = VideoWriter(avifilename);
            open(vObj);
            writeVideo(vObj, mov);
            close(vObj);
        else
            aviobj = avifile(avifilename, 'fps', fps); %#ok<DAVIFL>
            for framenum = 1:length(mov)
                aviobj = addframe(aviobj, mov(framenum));
            end
            aviobj = close(aviobj); %#ok<NASGU>
        end
    end
else
    % Plot the color grid, using the minimum and maximum mean grid values
    % to reset the color scale:
    cmin = Inf;
    cmax = -Inf;
    for xidx = 1:numx
        for yidx = 1:numy
            if ~isempty(colorgrid{xidx, yidx}) && ...
                    length(colorgrid{xidx, yidx}) >= N
                if sum(isnan(colorgrid{xidx, yidx})) / ...
                        numel(colorgrid{xidx, yidx}) > 0.5
                    warning( 'lfp_disp2d:NaNs', ...
                        'The grid bin starting at (x,y) = (%d, %d) contains %d/%d NaNs', ...
                        xvals(xidx), yvals(yidx), ...
                        sum(isnan(colorgrid{xidx, yidx})), ...
                        numel(colorgrid{xidx, yidx}) );
                end
                meanvalue = nanmean(colorgrid{xidx, yidx});
                cmin = min([cmin meanvalue]);
                cmax = max([cmax meanvalue]);
                patch( ...
                    [xvals(xidx) xvals(xidx) xvals(xidx+1) xvals(xidx+1)], ...
                    [yvals(yidx) yvals(yidx+1) yvals(yidx+1) yvals(yidx)], ...
                    meanvalue, 'Parent', hA,...
                    'CDataMapping', 'scaled', 'FaceColor', 'flat', ...
                    'EdgeColor', 'none' );
            end
        end
    end
    if isinf(cmin)
        % Then <cmax> is too, and conversely
        error('lfp_disp2d:nodata', ...
            'There are no ''colorgrid'' values to plot.');
    else
        caxis(hA, [cmin cmax]);
    end
end
