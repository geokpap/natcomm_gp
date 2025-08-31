function dg_gunCSV2mat(filename, varargin)
% Reads a CSV file and creates a set of *.mat files that can be read
% directly into lfp_lib (*.evtsav and 'dg_Nlx2Mat' formats).
%INPUTS
%OUTPUTS
% All output goes directly to new files in the same directory as
% <filename>.  They are in Matlab '-v7.3' format.  The filenames are as
% follows, where <basename> is the base name of <filename>:
%   <basename>.evtsav: events and trial params.
%       Trial Parameters: 1-3 are boolean (but saved as double), 4 is
%           natively double. 
%       1. Airpuff (low-probability punishment in the
%       incorrect port)
%       2. Left port is correct (if <false>, then Right port is correct)
%       3. Is switch trial
%       4. Value (calculated from learning model)
%   <basename>RedL.mat: red left-side photometry channel.
%   <basename>RedR.mat: red right-side photometry channel.
%   <basename>GreenL.mat: green left-side photometry channel.
%   <basename>GreenR.mat: green right-side photometry channel.
%OPTIONS
% 'dest', destdir - writes output to <destdir> instead of to the same
%   directory as the source file.
%NOTES
% Event columns can start anywhere in <filename>, but we assume that the
% last event column is "Right error", and that the immediately successive
% columns are "Airpuff", "Left or Right correct(Left=True)", "L/R switch
% count".  The column immediately before the first event column
% ("Entries(even if out of task)") is assumed to be the timestamp column.
% The column number for "value" is semi-hard-coded as <valcol>.  The
% photometry columns are semi-hard-coded as <photcols>.

%$Rev: 306 $
%$Date: 2023-10-05 16:32:32 -0400 (Thu, 05 Oct 2023) $
%$Author: dgibson $

[destdir, basename, ~] = fileparts(filename);

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dest'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        otherwise
            error('dg_gunCSV2mat:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Event ID codes:
InitOn = 1; % initiation port
InitOff = 2;
LOn = 3; % left port
LOff = 4;
ROn = 5; % right port
ROff = 6;
TrialStart = 7;
TrialEnd = 8;
TurnStart = 9;
TurnEnd = 10;
LRwd = 11;
RRwd = 12;
Initiation = 13;
LNoRwd = 14;
RNoRwd = 15;
LError = 16;
RError = 17;

secPerDay = 24 * 60 * 60;
allocsize = 2^12;
events = NaN(allocsize, 2);
tp = {};
evtidx = 0;
trialnum = 0;
chnames = {'RedL' 'RedR' 'GreenL' 'GreenR'};

cells = dg_tabread3(filename, 'delim', ',');
eventscols = find(ismember(cells(1, :), 'Entries(even if out of task)')) ...
    : find(ismember(cells(1, :), 'Right error'));
TScol = eventscols(1) - 1;
valcol = 68;
photcols = 57:60;
t0 = datenum(cells{2, TScol}); %#ok<DATNM> 

oldevents = [];
for linenum = 2:size(cells,1)
    currentevents = NaN(1, length(eventscols));
    for evtcolidx = [1:3, length(eventscols)]
        currentevents(evtcolidx) =  isequal( ...
            cells{linenum, eventscols(evtcolidx)}, 'True' );
    end
    currentevents(4:end-1) = str2double( ...
        cells(linenum, eventscols(4:end-1)) );
    if isempty(oldevents)
        oldevents = currentevents;
        continue
    end
    for chgidx = find(currentevents ~= oldevents)
        switch chgidx
            case 1 % 'Entries(even if out of task)'
                if currentevents(chgidx)
                    addevent(InitOn, cells{linenum, TScol});
                else
                    addevent(InitOff, cells{linenum, TScol});
                end
            case 2 % 'left entry'
                if currentevents(chgidx)
                    addevent(LOn, cells{linenum, TScol});
                else
                    addevent(LOff, cells{linenum, TScol});
                end
            case 3 % 'right entry'
                if currentevents(chgidx)
                    addevent(ROn, cells{linenum, TScol});
                else
                    addevent(ROff, cells{linenum, TScol});
                end
            case 4 % 'Trial On'
                if currentevents(chgidx)
                    addevent(TrialStart, cells{linenum, TScol});
                    trialnum = trialnum + 1;
                else
                    addevent(TrialEnd, cells{linenum, TScol});
                end
            case 5 % 'In Turn Area'
                if currentevents(chgidx)
                    addevent(TurnStart, cells{linenum, TScol});
                else
                    addevent(TurnEnd, cells{linenum, TScol});
                end
            case 6 % 'Left correct+rewarded'
                addevent(LRwd, cells{linenum, TScol});
            case 7 % 'Right correct+rewarded'
                addevent(RRwd, cells{linenum, TScol});
            case 8 % 'Initiation(start) counter'
                addevent(Initiation, cells{linenum, TScol});
                tp{1, trialnum}(1, 1) = str2double( ...
                    cells{linenum, eventscols(end) + 1} ); %#ok<*AGROW> 
                tp{1, trialnum}(1, 2) = str2double( ...
                    cells{linenum, eventscols(end) + 1} );
                tp{1, trialnum}(1, 3) = ~isequal( ...
                    cells{linenum, eventscols(end) + 1}, ...
                    cells{linenum - 1, eventscols(end) + 1} );
                % There may not be any <valcol>:
                if valcol > size(cells, 2)
                    tp{1, trialnum}(1, 4) = NaN;
                else
                    tp{1, trialnum}(1, 4) = str2double( ...
                        cells{linenum, valcol} );
                end
            case 9 % 'Left correct+omission'
                addevent(LNoRwd, cells{linenum, TScol});
            case 10 % 'Right correct+omission'
                addevent(RNoRwd, cells{linenum, TScol});
            case 11 % 'Left error'
                addevent(LError, cells{linenum, TScol});
            case 12 % 'Right error'
                addevent(RError, cells{linenum, TScol});
            otherwise
                error('dg_gunCSV2mat:chgidx', ...
                    '<chgidx> out of range at line %d', linenum);
        end
    end
    oldevents = currentevents;
end
lfp_save_events = events(~isnan(events(:,1)), :);
lfp_save_params = tp;
OutputFileName = sprintf('%s.evtsav', basename);
save(fullfile(destdir, OutputFileName), ...
    'lfp_save_events', 'lfp_save_params');
dg_Nlx2Mat_Timestamps = round(1e6 * ...
    (datenum(cells(2:end, TScol)) - t0) * secPerDay ); %#ok<DATNM> 
dg_Nlx2Mat_SamplesUnits = 'SDs';
for chidx = 1:length(photcols)
    dg_Nlx2Mat_Samples = reshape( ...
        str2double(cells(2:end, photcols(chidx))), ...
        1, [] );
    OutputFileName = sprintf('%s%s.mat', basename, chnames{chidx});
    save(fullfile(destdir, OutputFileName), ...
        'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_Samples', ...
        'dg_Nlx2Mat_SamplesUnits', '-v7.3' );
end

    function addevent(evtID, timestamp)
        %INPUTS
        % evtID: see "Event ID codes".
        % timestamp: string containing date-time code.
        evtidx = evtidx + 1;
        if evtidx > size(events, 1)
            events = [ events
                NaN(allocsize, 2) ];
        end
        events(evtidx, :) = [(datenum(timestamp) - t0) * secPerDay, evtID]; %#ok<DATNM> 
    end

end

