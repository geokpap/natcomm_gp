function hF = dg_viewTfile(filename, clustnum, varargin)
%dg_viewTfile(filename, clustnum)
% Displays a simple raster plot showing spikes from clustnum with events
% read from the T file <filename>.  <hF> is a handle to the newly created
% figure.  You can click on the event markers to get event details
% displayed on the command window.
%OPTIONS
% 'allowedpairs', allowedpairs - a list of exceptions to the rule that
%   there must only be one match to each element of <eventlist> (see
%   below).  Does nothing without 'eventlist'.  <allowedpairs> is a cell
%   array, each of whose elements contains a list of event IDs that are
%   allowed to occur within the same trial even though they belong to the
%   same element of <eventlist>.  The lists can actually be longer than
%   just pairs, and as long as the events that match a single cell of
%   <eventlist> are all listed together in the same cell of <allowedpairs>,
%   no error is raised.
% 'eventlist', eventlist - <eventlist> is a cell array where each entry
%   is a row vector of possible alternative event IDs that should exist in
%   each trial; if one is missing a warning is raised. They must also occur
%   in the same order in which they are listed in <eventlist>; otherwise,
%   an error is raised after constructing the plot. An error is also raised
%   if there is more than one event that matches any single cell in
%   <eventlist>, unless the pair of matched events is listed on
%   <allowedpairs> (see above). 
% 'evts2skip', evts2skip - <evts2skip> is a list of event IDs that should
%   be totally ignored.  It defaults to 50, which normally contains
%   garbage.
% 'maxdur', maxdur - if any event or spike in the T file has a timestamp
%   greater than <maxdur>, an error is raised.
% 'noplot' - too complicated to understand or even to describe.
% 'warningsonly' - raises no recoverable errors, but only fatal errors and
%   warnings 

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

allowedpairs = [];
eventlist = [];
evts2skip = 50;
maxdur = 0;
plotflag = true;
warningsonlyflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'allowedpairs'
            argnum = argnum + 1;
            allowedpairs = varargin{argnum};
        case 'eventlist'
            argnum = argnum + 1;
            eventlist = varargin{argnum};
        case 'evts2skip'
            argnum = argnum + 1;
            evts2skip = varargin{argnum};
        case 'maxdur'
            argnum = argnum + 1;
            maxdur = varargin{argnum};
        case 'noplot'
            plotflag = false;
        case 'warningsonly'
            warningsonlyflag = true;
        otherwise
            error('dg_viewTfile:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

colors = {'r' 
    'g' 
    'b' 
    'm'
    'c'
    [.7 .7 0]
    [0 .7 .7]
    [.7 0 .7]
    [.5 .5 .5]
    };

[fileheader, trialdata] = ...
    dg_ReadRodentFormat(filename);
if fileheader.TSize == 0
    if warningsonlyflag
        warning('dg_viewTfile:notrials', ...
            'The T file %s contains zero trials.', filename);
    else
        error('dg_viewTfile:notrials', ...
            'The T file %s contains zero trials.', filename);
    end
end
if clustnum > fileheader.CSize
    error('dg_viewTfile:clustnum', ...
        '%s : Cannot show cluster %d - there are only %d clusters in the file.', ...
        filename,clustnum, fileheader.CSize);
end
eventses = cell(size(trialdata));
[eventses{:}] = deal(trialdata.events);
eventses = cell2mat(eventses');
eventses(:, evts2skip) = 0;
latespikes = cell(size(trialdata));

if plotflag
    hF = figure;
    hA = axes;
    hold(hA, 'on');
    axis ij;
    grid on;
    ylabel('Trial Number');
    xlabel('Clock ticks (0.1 ms)');
    title(sprintf('%s cluster %d', filename, clustnum), ...
        'Interpreter', 'none');
end

for trial = 1:length(trialdata)
    spikes = trialdata(trial).spikes(:,clustnum);
    if plotflag
        dg_plottick(spikes(spikes>0), trial, 1, 'k');
    end
    if maxdur
        latespikes{trial} = find(spikes>maxdur);
    end
    for evtix = 1:size(eventses,2)
        if eventses(trial, evtix)
            if plotflag
                eventcolor = colors{mod(evtix, length(colors)) + 1};
                hL = plot(hA, ...
                    eventses(trial, evtix), ...
                    trial, ...
                    'Color', eventcolor, ...
                    'Marker', '.' );
                if isequal(class(eventcolor), 'char')
                    eventcolorstr = eventcolor;
                else
                    eventcolorstr = mat2str(eventcolor);
                end
                detailstr = sprintf( ...
                    '\\nTrial=%d\\nTimestamp=%.0f\\nColor="%s"\\nEventID=%.0f (0x%X)', ...
                    trial, eventses(trial, evtix), ...
                    eventcolorstr, evtix, evtix );
                set(hL, ...
                    'ButtonDownFcn', ...
                    ['fprintf(1,''' detailstr '\n'')'] );
            end
        end
    end
end

for trial = 1:length(trialdata)
    evtsintrial = find(eventses(trial,:));
    evtIDs = cell(size(eventlist));
    for k = 1:length(eventlist)
        if ~any(ismember(evtsintrial, eventlist{k}))
            warning('dg_viewTfile:evtmissing', ...
                '%s: Trial %d has no event %s', ...
                filename, trial, mat2str(eventlist{k}) );
        end
        evtIDs{k} = evtsintrial(ismember(evtsintrial, eventlist{k}));
        if length(evtIDs{k}) > 1
            OK = false;
            for allowedidx = 1:length(allowedpairs)
                if all(ismember(evtIDs{k}, allowedpairs{allowedidx}))
                    OK = true;
                    break
                end
            end
            if ~OK
                if warningsonlyflag
                    warning('dg_viewTfile:eventdup', ...
                        '%s: Trial %d has redundant events %s', ...
                        filename, trial, mat2str(evtIDs{k}) );
                else
                    error('dg_viewTfile:eventdup', ...
                        '%s: Trial %d has redundant events %s', ...
                        filename, trial, mat2str(evtIDs{k}) );
                end
            end
        end
        if k > 1 
            %&& eventses(trial, evtIDs{k}) <= eventses(trial, evtID(k-1))
            for evtIDidx = 1:length(evtIDs{k})
                for prevevtIDidx = 1:length(evtIDs{k-1})
                    if eventses(trial, evtIDs{k}(evtIDidx)) <= ...
                            eventses(trial, evtIDs{k-1}(prevevtIDidx))
                        if warningsonlyflag
                            warning('dg_viewTfile:evtdisordered', ...
                                '%s: Trial %d has event %d before (or same time as) event %d', ...
                                filename, trial, evtIDs{k}(evtIDidx), evtIDs{k-1}(prevevtIDidx) );
                        else
                            error('dg_viewTfile:evtdisordered', ...
                                '%s: Trial %d has event %d before (or same time as) event %d', ...
                                filename, trial, evtIDs{k}(evtIDidx), evtIDs{k-1}(prevevtIDidx) );
                        end
                    end
                end
            end
        end
    end
    if maxdur && any(any(eventses > maxdur))
        msg = 'The following trials have events > maxdur:';
        for idx = reshape(find(eventses > maxdur), 1, [])
            [I,J] = ind2sub(size(eventses), idx);
            msg = sprintf('%s\nTrial %d evt %d: %.0f', msg, ...
                I, J, eventses(I,J) );
        end
        if warningsonlyflag
            warning('dg_viewTfile:maxdur1', msg);
        else
            error('dg_viewTfile:maxdur1', msg);
        end
    end
    if ~isempty(latespikes{trial})
        if warningsonlyflag
            warning('dg_viewTfile:maxdur2', ...
                'Trial %d has spikes %s > maxdur in time range [%.0f %.0f]', ...
                trial, dg_canonicalSeries(latespikes{trial}), ...
                trialdata(trial).spikes(latespikes{trial}(1),clustnum), ...
                trialdata(trial).spikes(latespikes{trial}(end),clustnum) );
        else
            error('dg_viewTfile:maxdur2', ...
                'Trial %d has spikes %s > maxdur in time range [%.0f %.0f]', ...
                trial, dg_canonicalSeries(latespikes{trial}), ...
                trialdata(trial).spikes(latespikes{trial}(1),clustnum), ...
                trialdata(trial).spikes(latespikes{trial}(end),clustnum) );
        end
    end
end
