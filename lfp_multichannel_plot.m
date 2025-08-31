function [hF, data] = lfp_multichannel_plot(plothandle, figtitle, plots, ...
    sampledata, auxdata, reftime, axesinfo, option, varargin)
%lfp_multichannel_plot(plothandle, figtitle, plots, ...
%    sampledata, auxdata, reftime, axesinfo, option, bigevts)
%lfp_multichannel_plot(plothandle, figtitle, plots, ...
%    sampledata, auxdata, reftime, axesinfo, option, 'dB')

%AKA the clever design concept that should never have been
%   DO NOT MODIFY THIS!  REPLACE IT INSTEAD!

% INPUTS
%  plothandle: the figure or axes handle into which to plot.  To create a
%       new figure, specify a new figure number, or specify [] to
%       let Matlab choose a figure number.  For some reason it was
%       originally necessary to destroy the old figure and create a new one
%       if <plothandle> was actually a figure handle rather than an
%       integer; however, in the case where 'Visible' is 'off', that
%       behavior is now overridden.
%  figtitle: title string for top plot; backslashes will automatically
%       be escaped so as not to confuse the TeX interpreter.
%  plots: In Standard mode (see sampledata), file numbers to plot; these 
%       should be  members of lfp_ActiveFilenums. In Custom mode, cell row
%       vector of strings for labelling each plot; must have same number of
%       columns as sampledata.
%  sampledata: 
%       If this contains a single row, then lfp_multichannel_plot runs in
%       Standard mode; if it contains multiple rows, then we run in Custom
%       mode. In Standard mode, sampledata contains the indices into
%       lfp_Samples{filenum} of the sample points to plot.  In Custom mode,
%       sampledata contains a matrix containing column vectors of literal
%       data to plot for one or more channels, one channel per column.  The
%       first row is considered to be absolute time 0, and absolute time
%       for subsequent rows is ((row-1) * lfp_SamplePeriod).  If there is a
%       third dimension to the array, then sampledata(:,:,2) is displayed
%       as error bars (except see notes on <option> below).
%  auxdata: 
%       In Standard mode: the indices in lfp_Events of the events
%           to plot.  
%       In Custom mode: arbitrary events in the format of lfp_Events,
%           i.e. [ timestamp, event ID], one row per event. Timestamps are
%           in same time frame as sampledata. Alternatively, auxdata may be
%           a column vector containing arbitrary values to use for the
%           x-coordinates of the plot; in this case, length(auxdata) must
%           be the same as size(sampledata,1). Or, auxdata can be [], in
%           which case it is ignored.
%  reftime: absolute time to use as time 0 in the plots, or more generally
%       in Custom mode, a number that is subtracted from all the values
%       along the x axis before using them in plotting.
%  axesinfo:  structure argument with fields 'xlabel', 'xlim',
%       'ylim'.  axesinfo.xlim and axesinfo.ylim are used to set the limits
%       on the axes (empty matrix specifies auto scaling); xlabel field is
%       used to label the x axis.  If axesinfo.ylim has multiple rows, then
%       it is assumed that the rows contain a complete set of Y limits for
%       every channel being displayed.  In this case, the value [0 0] is
%       used to specify automatic scaling.  Optional fields:
%           'trialnums' - contains the trialnums that are to be used
%       as a basis for providing ButtonDownFcn info.  If the length of
%       axesinfo.trialnums does not match the number of traces plotted,
%       then it is ignored.  
%           'trigtimes' - used together with 'trialnums' for ButtonDownFcn 
%       info; if it is absent or its length does not match, NaN is used for
%       the trigger time.
%           'plotflag' - default is true; if false, skips figure creation.
%  option:  if optional parameter <option> is 'errorcurves', then
%       sampledata(:,:,2) is displayed as continuous curves above and below
%       sampledata(:,:,1) instead of error bars; if 'asymmetric',
%       sampledata(:,:,2:3) contains two different curves for upper and
%       lower bounds; if 'ovr', the third dimension of sampledata is
%       interpreted as multiple plots to superimpose on a single set of
%       axes; if 'ovrchan' then instead of showing one subplot for each
%       channel (i.e. each column of sampledata), the channels are overlaid
%       on a single plot and a legend is used instead of a ylabel.  If
%       <option> is '', the effect is the same as not giving <option> at
%       all.  NOTE: as of 8-Dec-2007, 'ovr' was used exclusively to show
%       overlaid trials; however, since fewer than all of the enabled
%       trials might be displayed, I had to make the ButtonDownFcn code
%       test for that possibility.
%  bigevts:  if given and not empty, this invokes behavior equivalent
%       to the 'bigevts' option of lfp_spec.
% OUTPUTS
%  hF: handle to the newly created figure
%  data: if <sampledata> has fewer than 3 dimensions, a cell array
%       containing a 2-column array of X and Y coordinates for each graph
%       plotted; if option is 'errorcurves', includes a 3rd column for the
%       error; for 'asymmetric', includes a 4th column for the red curve.

%$Rev: 390 $
%$Date: 2017-06-21 17:37:26 -0400 (Wed, 21 Jun 2017) $
%$Author: dgibson $

global lfp_SamplePeriod lfp_Samples lfp_XTicksOnAll lfp_FileNames ...
    lfp_SamplesUnits lfp_Events lfp_SelectedEventIDs lfp_EventColors ...
    lfp_EventNames lfp_EventDefaultColor

dBflag = false;
data = {};
if ~isempty(plothandle)
    if verLessThan('matlab','8.4.0')
        plot_int = plothandle;
    else
        plot_int = get(plothandle, 'Number');
    end
end

% literals representing different modes of operation
STD = 1;
CUST = 2;

if nargin < 9
    bigevts = {};
end
argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'dB'
                dBflag = true;
                argnum = argnum + 1;
                dBref = varargin{argnum};
        end
    else
        bigevts = varargin{argnum};
        if ~isempty(bigevts)
            bigevtIDs = cell2mat(bigevts(:,1));
        end
    end
    argnum = argnum + 1;
end

if nargin < 7
    error('lfp_multichannel_plot:noaxesinfo', ...
        '<axesinfo> must be supplied.');
else
    for fname = {'xlabel' 'xlim' 'ylim'}
        if ~isfield(axesinfo, fname)
            error('lfp_multichannel_plot:badaxesinfo', ...
                [ 'The axesinfo argument must contain these fields: ' ...
                'xlabel, xlim, ylim.' ]);
        end
    end
end
if ~isfield(axesinfo, 'plotflag')
    axesinfo.plotflag = true;
end
if ~isempty(axesinfo.ylim) && (size(axesinfo.ylim, 2) ~= 2)
    warning('lfp_multichannel_plot:axesinfo.ylim', ...
        'axesinfo.ylim must be an Nx2 array of double.');
end

if nargin < 8
    option = '';
end

if size(sampledata,1) > 1
    mode = CUST;
    if length(plots) ~= size(sampledata,2)
        error('lfp_multichannel_plot:channelcountmismatch', ...
            'In Custom mode, plots must have the same number of columns as sampledata.');
    end
else
    mode = STD;
end

switch mode
    case STD
        if length(sampledata) < 2
            error('lfp_multichannel_plot:sampledata', ...
                'Only %d sample selected, not enough to plot', ...
                length(sampledata) );
        end
        if size(plots,1) ~= 1 ...
                || ~(strcmp(class(plots), 'double')) ...
                || ~all(fix(plots) == plots)
            error('lfp_multichannel_plot:badplots1', ...
                '<plots> must be an integer row vector.');
        end
    case CUST
        if size(sampledata, 1) < 2
            error('lfp_multichannel_plot:sampledata2', ...
                'Only %d sample selected, not enough to plot', ...
                length(sampledata) );
        end
        if size(plots,1) ~= 1 ...
                || ~(strcmp(class(plots), 'cell')) ...
                || ~(strcmp(class(plots{1}), 'char'))
            error('lfp_multichannel_plot:badplots2', ...
                '<plots> must be a cell vector of strings.');
        end
end

if ~isempty(plothandle)
    if ishandle(plothandle)
        if isequal(get(plothandle, 'Type'), 'figure')
            if isequal(get(plothandle, 'Visible'), 'off')
                % This figure handle was specially created by
                % the caller.
                hF = plothandle;
            else
                % close any previous incarnation of the figure:
                close(plothandle);
                plothandle = figure(plot_int);
                hF = plothandle;
            end
        else
            if ~isequal(get(plothandle, 'Type'), 'axes')
                error('lfp_multichannel_plot:plothandle', ...
                    '<plothandle> must be a handle to a figure or axes.' );
            end
            hF = get(plothandle, 'Parent');
        end
    else
        hF = figure(plothandle);
    end
elseif axesinfo.plotflag
    hF = figure;
    plothandle = hF;
else
    hF = [];
end
plotnum = 0;
switch mode
    case STD
        xscale = lfp_index2time(sampledata) - reftime;
        channellist = plots;
        auxdata = reshape(auxdata, 1, []);
    case CUST
        if isempty(auxdata) || size(auxdata,2) == 2
            xscale = (0 : lfp_SamplePeriod : ...
                ((size(sampledata,1) - 1) * lfp_SamplePeriod) ) - reftime;
        elseif size(auxdata,2) == 1
            if length(auxdata) ~= size(sampledata,1)
                error('lfp_multichannel_plot:badxscale', ...
                    'length(auxdata) does not match sampledata');
            end
            xscale = auxdata';
        else
            error('lfp_multichannel_plot:badauxdata', ...
                'auxdata has an illegal number of columns');
        end
        channellist = 1:size(sampledata,2);
end
if isequal(option, 'ovrchan')
    numplots = 1;
    channellist = 1;
else
    numplots = length(plots);
end
for channel = channellist
    plotnum = plotnum + 1;
    if ~isempty(plothandle) && ishandle(plothandle) && isequal(get(plothandle, 'Type'), 'figure')
        subplot(numplots, 1, plotnum);
    end
    if ~isempty(hF)
        hold on;
    end
    
    % Plot the CSC data:
    hL = [];
    switch mode
        case STD
            if isequal(option, 'ovrchan')
                samples = zeros(numel(sampledata), numel(plots));
                for k = 1:numel(plots)
                    samples(:,k) = lfp_Samples{plots(k)}(sampledata);
                end
                if axesinfo.plotflag
                    hL = plot(xscale, samples);
                end
                for k=1:length(hL)
                    set(hL(k), 'ButtonDownFcn', ...
                        sprintf('disp(''trace %d'')', k));
                end
            else
                samples = lfp_Samples{channel};
                data{plotnum} = ...
                    [xscale' reshape(samples(sampledata), [], 1)];
                if axesinfo.plotflag
                    plot(xscale, samples(sampledata));
                end
            end
            if axesinfo.plotflag && size(axesinfo.ylim,1) > 1
                if any(axesinfo.ylim(channel,:))
                    ylim(axesinfo.ylim(channel,:));
                end
            end
        case CUST
            if size(sampledata, 3) > 1 || isequal(option, 'ovrchan')
                % Multiple traces per plot
                switch option
                    case 'ovr'
                        if axesinfo.plotflag
                            hL = plot(xscale, ...
                                squeeze(sampledata(:,channel,:))');
                        end
                        % If the number of traces is not equal to the
                        % number of enabled trials, then we still label
                        % them as 'trace %d'.
                        if isfield(axesinfo, 'trialnums')
                            trials2plot = axesinfo.trialnums;
                        else
                            trials2plot = lfp_enabledTrials;
                        end
                        istrials = length(hL) == length(trials2plot);
                        if istrials
                            if ~isfield(axesinfo, 'trigtimes') || ...
                                    length(axesinfo.trigtimes) ~= ...
                                    length(trials2plot)
                                trigtimes = NaN(size(trials2plot));
                            else
                                trigtimes = axesinfo.trigtimes;
                            end
                            for k=1:length(hL)
                                set(hL(k), 'ButtonDownFcn', ...
                                    sprintf('disp(''%s trigtime=%.6f''); lfp_declareGlobals; lfp_ClickedTrials = union(lfp_ClickedTrials, %d);', ...
                                    lfp_getTrialID(trials2plot(k)), ...
                                    trigtimes(k), ...
                                    trials2plot(k) ));
                            end
                        else
                            for k=1:length(hL)
                                set(hL(k), 'ButtonDownFcn', ...
                                    sprintf('disp(''trace %d'')', k));
                            end
                        end
                    case 'ovrchan'
                        if axesinfo.plotflag
                            hL = plot(xscale, squeeze(sampledata(:,:,1))');
                        end
                        for k=1:length(hL)
                            set(hL(k), 'ButtonDownFcn', ...
                                sprintf('disp(''trace %d'')', k));
                        end
                    case 'errorcurves'
                        data{plotnum} = [ ...
                            reshape(xscale, [], 1) ...
                            reshape(sampledata(:,channel,1), [], 1) ...
                            reshape(sampledata(:,channel,2), [], 1)];
                        if axesinfo.plotflag
                            if dBflag
                                plot(xscale, ...
                                    10*log10(sampledata(:,channel,1)') - dBref, 'k');
                                plot(xscale, ...
                                    10*log10(sampledata(:,channel,1)' + ...
                                    sampledata(:,channel,2)') - dBref, 'r');
                                plot(xscale, ...
                                    10*log10(sampledata(:,channel,1)' - ...
                                    sampledata(:,channel,2)') - dBref, 'b');
                            else
                                plot(xscale, sampledata(:,channel,1)', 'k');
                                plot(xscale, sampledata(:,channel,1)' + ...
                                    sampledata(:,channel,2)', 'r');
                                plot(xscale, sampledata(:,channel,1)' - ...
                                    sampledata(:,channel,2)', 'b');
                            end
                        end
                    case 'asymmetric'
                        data{plotnum} = [ ...
                            reshape(xscale, [], 1) ...
                            reshape(sampledata(:,channel,1), [], 1) ...
                            reshape(sampledata(:,channel,2), [], 1) ...
                            reshape(sampledata(:,channel,3), [], 1)];
                        if axesinfo.plotflag
                            plot(xscale, sampledata(:,channel,1)', 'k');
                            plot(xscale, sampledata(:,channel,2)', 'b');
                            plot(xscale, sampledata(:,channel,3)', 'r');
                        end
                    otherwise
                        errorbar(xscale, sampledata(:,channel,1)', ...
                            sampledata(:,channel,2)');
                end
            else
                % Single trace per plot
                data{plotnum} = [xscale' sampledata(:,channel,1)];
                if axesinfo.plotflag
                    plot(xscale, sampledata(:,channel,1)');
                end
            end
            if axesinfo.plotflag && size(axesinfo.ylim,1) > 1
                if any(axesinfo.ylim(channel,:))
                    ylim(axesinfo.ylim(channel,:));
                end
            end
    end
    
    if axesinfo.plotflag
        hA = gca;
        if plotnum == 1
            title(figtitle, 'Interpreter', 'none');
        end
        if plotnum == numplots
            xlabel(axesinfo.xlabel);
        else
            if ~lfp_XTicksOnAll
                set(hA, 'XTickLabel', '');
            end
        end
        grid on;
        switch mode
            case STD
                if isequal(option, 'ovrchan')
                    hleg = legend(hL, lfp_FileNames{plots});
                    set(hleg, 'Interpreter', 'none');
                else
                    ylabel(sprintf('%s, %s', lfp_FileNames{channel}, ...
                        lfp_SamplesUnits{channel}), 'Interpreter', 'none');
                end
                eventlist = auxdata;
                events = lfp_Events;
            case CUST
                if isequal(option, 'ovrchan')
                    hleg = legend(hL, plots{:});
                    set(hleg, 'Interpreter', 'none');
                else
                    ylabel(plots{channel}, 'Interpreter', 'none');
                end
                eventlist = 1:size(auxdata,1);
                events = auxdata;
        end
        if mode == STD || (mode == CUST && size(auxdata,2) == 2)
            % Plot the event data:
            for evtix = eventlist
                eventtime = events(evtix, 1) - reftime;
                evtid = events(evtix,2);
                if evtid == 0
                    warning('lfp_multichannel_plot:eventID0', ...
                        'Event at %d has ID 0', eventtime );
                else
                    if evtid == 0
                        warning('lfp_multichannel_plot:eventID0', ...
                            'Event number %d has ID 0', evtix );
                    elseif (evtid <= length(lfp_SelectedEventIDs)) ...
                            && lfp_SelectedEventIDs(evtid) ...
                            || (evtid > length(lfp_EventColors)) ...
                            && ~isempty(bigevts) ...
                            && ismember(evtid, bigevtIDs)
                        % Either the event ID is selected, or it is a
                        % specified big event ID.
                        eventcolor = '';
                        eventname = ''; % required for detailstr
                        if evtid <= length(lfp_EventNames)
                            eventname = lfp_EventNames{evtid};
                        end
                        if evtid <= length(lfp_EventColors)
                            eventcolor = lfp_EventColors{evtid};
                        else
                            for k = 1:size(bigevts,1)
                                if ismember(evtid, bigevts{k,1}')
                                    eventcolor = bigevts{k,2};
                                    break
                                end
                            end
                        end
                        if isempty(eventcolor)
                            eventcolor = lfp_EventDefaultColor;
                        end
                        if axesinfo.plotflag
                            hL = plot([ eventtime eventtime ], ...
                                get(hA, 'YLim'), ...
                                'Color', eventcolor );
                        end
                        if mode == STD
                            eventTS = lfp_Events(evtix,1);
                        else
                            eventTS = eventtime;
                        end
                        if exist('eventcolor', 'var') ...
                                && exist('eventname', 'var') ...
                                && exist('evtid', 'var')
                            if isequal(class(eventcolor), 'char')
                                eventcolorstr = eventcolor;
                            else
                                eventcolorstr = mat2str(eventcolor);
                            end
                            detailstr = sprintf( ...
                                '\\nTimestamp=%.6f\\nLineColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
                                eventTS, ...
                                eventcolorstr, eventname, evtid, evtid );
                        else
                            detailstr = 'No details available';
                        end
                        set(hL, ...
                            'ButtonDownFcn', ...
                            ['fprintf(1,''' detailstr '\n'')'] );
                    end
                end
            end
        end
    end
end

% Apply axesinfo.xlim:
if axesinfo.plotflag && ~isempty(axesinfo.xlim)
    if ~all(size(axesinfo.xlim) == [1 2]) ...
            || ~strcmp(class(axesinfo.xlim), 'double')
        warning('lfp_multichannel_plot:badxlimall', ...
            'axesinfo.xlim must be a 1x2 array of double.');
    else
        dg_xlimall(plothandle, axesinfo.xlim);
        xlimidx = find(xscale >= axesinfo.xlim(1) ...
            & xscale <= axesinfo.xlim(2) );
        for plotnum = 1:length(data)
            data{plotnum} = data{plotnum}(xlimidx, :);
        end
    end
end

% Apply axesinfo.ylim to all plots:
if axesinfo.plotflag && ~isempty(axesinfo.ylim) && ...
        isequal(size(axesinfo.ylim), [1 2]) ...
        && axesinfo.ylim(2) > axesinfo.ylim(1)
    dg_ylimall(plothandle, axesinfo.ylim);
end
if axesinfo.plotflag
    hold off;
end


