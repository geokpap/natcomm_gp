function [xposn, yposn, outliersx, outliersy, medians] = ...
    lfp_getEvtPosns(x, y, evtIDs, varargin)
%OPTIONS
% 'trials', trials - <trials> works as usual

%$Rev: 136 $
%$Date: 2010-07-02 23:08:22 -0400 (Fri, 02 Jul 2010) $
%$Author: dgibson $

lfp_declareGlobals;

trials = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'trials'
            argnum = argnum + 1;
            trials = varargin{argnum};
        otherwise
            error('lfp_getEvtPosns:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if isempty(trials)
    trials = lfp_enabledTrials(find(lfp_SelectedTrials));
else
    trials = lfp_enabledTrials(trials);
end

measuredevts = NaN(length(evtIDs), length(trials));
% Collect relevant event indices from trials(trialidx) periods
for trialidx = 1:length(trials)
    trialevents = lfp_Events( lfp_TrialIndex(trials(trialidx),1) ...
        : lfp_TrialIndex(trials(trialidx),2), : ); %#ok<NODEF>
    for evtIDidx = 1:length(evtIDs)
        evtidx = find(ismember(trialevents(:,2), evtIDs{evtIDidx}));
        if ~isempty(evtidx)
            measuredevts(evtIDidx, trialidx) = evtidx(1) + ...
                lfp_TrialIndex(trials(trialidx),1) - 1;
        end
    end
end

measuredevtsidx = find(~isnan(measuredevts));

% <measuredevtsamples>, like <measuredevts>, is in evtID x trial form.
measuredevtsamples = NaN(length(evtIDs), length(trials));
measuredevtsamples(measuredevtsidx) = lfp_time2index( ...
    lfp_Events(measuredevts(measuredevtsidx), 1) );
xposn = NaN(length(evtIDs), length(trials));
yposn = NaN(length(evtIDs), length(trials));
for trialidx = 1:length(trials)
    % <xdata> and <ydata> have one row for each event in a single trial,
    % and one column for each sample in the vicinity of the event.
    samples = repmat((-2:2), length(evtIDs), 1) + ...
        repmat(measuredevtsamples(:, trialidx), 1, 5);
    isnotnan = ~isnan(samples);
    xdata = NaN(size(samples));
    xdata(isnotnan) = lfp_Samples{x}(samples(isnotnan)); %#ok<USENS>
    ydata = NaN(size(samples));
    ydata(isnotnan) = lfp_Samples{y}(samples(isnotnan));
    % xposn and yposn are deliberately NOT calculated nan-tolerantly; i.e.
    % we insist on having all five points or we don't include this trial
    % for that event. Like <measuredevtsamples> and <measuredevts>, they
    % are evtID x trial.
    xposn(:, trialidx) = mean(xdata, 2);
    yposn(:, trialidx) = mean(ydata, 2);
end
% Throw out data where either of the coordinates is NaN:
trialeventisnan = isnan(xposn) | isnan(yposn);
xposn(trialeventisnan) = NaN;
yposn(trialeventisnan) = NaN;

% For each event ID, keep removing outliers that are more than 20 pixels
% from the median until there are none left.  Maintain list of outlier
% coordinates to plot.
outliersx = NaN(size(xposn));
outliersy = NaN(size(yposn));
for evtIDidx = 1:length(evtIDs)
    tryagain = true;
    while tryagain
        % medians is in evtID x axis form, where col 1 is x and
        % col 2 is y. 
        medians(evtIDidx, 1) = ...
            dg_nanTolerantMedian(xposn(evtIDidx, :), 2);
        medians(evtIDidx, 2) = ...
            dg_nanTolerantMedian(yposn(evtIDidx, :), 2);
        if any(isnan(medians(evtIDidx, :)))
            % Give up immediately
            tryagain = false;
            warning('lfp_measureTMaze2:badoutliers', ...
                'Could not get rid of outliers for event %s', ...
                mat2str(evtIDs{evtIDidx}) );
        else
            xbad = abs(xposn(evtIDidx,:) - medians(evtIDidx,1)) > 20;
            ybad = abs(yposn(evtIDidx,:) - medians(evtIDidx,2)) > 20;
            isoutlier = xbad | ybad;
            if any(isoutlier)
                msg = sprintf( ...
                    'Throwing out outliers for event %s, trials %s', ...
                    mat2str(evtIDs{evtIDidx}), mat2str(trials(isoutlier)) );
                lfp_log(['lfp_measureTMaze2: ' msg]);
                warning('lfp_measureTMaze2:outliers', '%s', msg);
                outliersx(evtIDidx,isoutlier) = xposn(evtIDidx,isoutlier);
                outliersy(evtIDidx,isoutlier) = yposn(evtIDidx,isoutlier);
                xposn(evtIDidx,isoutlier) = NaN;
                yposn(evtIDidx,isoutlier) = NaN;
            else
                % No more outliers
                tryagain = false;
            end
        end
    end
end