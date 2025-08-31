function [table, header] = lfp_theresaEventTabulation(mode)
%LFP_THERESAEVENTTABULATION produces a table summarizing trial events
%[table, header] = lfp_theresaEventTabulation(mode)
%lfp_theresaEventTabulation(mode)

%[table, header] = lfp_theresaEventTabulation(mode)
%  Returns a table summarizing events for one task type, with one row for
%  each enabled trial. The task type is controlled by <mode>, which may be:
%   'fixation'
%   'scan'
%   'all'
%  <table> is a cell column vector, where each cell contains a numerical
%  row vector of variable length. <header> is a cell string row vector
%  containing headers for the columns of the elements in <table>.  The
%  "trial" column contains internal trial numbers.
%  NOTE: 'RexnT' or (for 'fixation' trials) 'FixDur' may be given as 0,
%  indicating that the required data were missing from the trial.  Note
%  that any events that take place between lfp_NominalTrialEnd and the
%  following lfp_NominalTrialStart are not included.
%
%lfp_theresaEventTabulation(mode)
%  Creates a tab-delimited text file containing the same table, with
%  headers on the first row.  The "trial" column contains Unique Trial IDs
%  instead of internal trial numbers.

% NOTE: this code assumes that the lfp_TrialParams are unmodified, i.e.
% that the high order byte of each one is the Param number, and the low
% order byte is the value.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

fixflag = false;
scanflag = false;
table = {};

lfp_getEvtIDs_theresa;
trial = 1;  % solely to keep lfp_getBDParams_theresa happy.
lfp_getBDParams_theresa;

% Initialize eventnum to spotnum table
numTargetRows = 5;
for spotnum = 1:numTargetRows^2
    pos(hex2dec('20') + spotnum) = spotnum;
end
for spotnum = 1:numTargetRows^2
    BDpos(spotnum) = spotnum;
end

% Initialize spotnum to [x, y] coordinates table - coords are Cartesian,
% i.e. x increases to the right, y increases upwards, and the lower left
% corner spot is located at (1, 1).
spotnums = flipud(reshape(1:numTargetRows^2, numTargetRows, numTargetRows));
for row = 1:numTargetRows
    for col = 1:numTargetRows
        coords(spotnums(row,col),:) = [col row];
    end
end

% Initialize headers
switch mode
    case 'fixation'
        fixflag = true;
        header = {'Trial', 'RexnT', 'Type', 'Result', ...
                'Subtype', 'FixPos', 'FixDur'};
    case 'scan'
        scanflag = true;
        header = {'Trial', 'RexnT', 'Type', 'Result', ...
                'RwdPos', 'ScanDur',...
                'Pos', 'Dur', 'OffDur', 'deltaX', 'deltaY', 'Pos', 'Dur'};
    case 'all'
        fixflag = true;
        scanflag = true;
        header = {'Trial', 'RexnT', 'Type', 'Result'};
    otherwise
        error('lfp_theresaEventTabulation:badmode', ...
            'Unknown mode "%s".', mode);
end

if nargout == 0
    OutputFileName = ['b_' mode '.xls'];
    [OutputFileName, OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, OutputFileName) );
    if isequal(OutputFileName, 0)
        return
    end
    outfid = fopen(fullfile(OutputPathName, OutputFileName), 'w');
    fprintf(outfid, '%s', header{1});
    fprintf(outfid, '\t%s', header{2:end});
    fprintf(outfid, '\n');
end

% Construct the table, which contains all the common (across tasks)
% parameters first, then task-specific parameters
for trial = lfp_enabledTrials(1:size(lfp_TrialIndex,1))
    datagood = true;
    trialevents = lfp_Events(lfp_TrialIndex(trial, 1) : ...
        lfp_TrialIndex(trial,2), :);
    trialtype = bitand(255, lfp_TrialParams{trial}(trialtypePN));
    trialresult = bitand(255, lfp_TrialParams{trial}(trialresultPN));
    if fixflag && (trialtype ==  fixationTrial)
        rexnt = trialevents(find(trialevents(:,2) == fixationOnsetID), 1) ...
            - trialevents(find(trialevents(:,2) == fixationStimOnID), 1);
        if isempty(rexnt)
            rexnt = 0;
        end
        duration = ...
            trialevents(find(trialevents(:,2) == fixationSpotOffID), 1) ...
            - trialevents(find(trialevents(:,2) == fixationOnsetID), 1);
        if isempty(duration)
            duration = 0;
        end
        % 'Subtype', 'FixPos', 'FixDur'
        taskspecific = [ bitand(255, lfp_TrialParams{trial}(subtypePN)) ...
            BDpos(bitand(255, lfp_TrialParams{trial}(fixPosPN))) ...
            duration ];
    elseif scanflag && (trialtype == scanTrial)
        % Get rid of extraneous events outside the period of interest
        scanAllTargetsOnindex = find(trialevents(:,2) == scanAllTargetsOnID);
        scanAllTargetsOffindex = find(trialevents(:,2) == scanAllTargetsOffID);
        if scanAllTargetsOnindex < length(trialevents)
            trialevents(1 : scanAllTargetsOnindex-1, :) = [];
        end
        if scanAllTargetsOffindex < length(trialevents)
            trialevents(scanAllTargetsOffindex+1 : end, :) = [];
        end
        % These are all column vectors (or matrices with time down
        % dimension 1)
        fixationstarts = find(ismember(trialevents(:, 2), scanTargetEnterID));
        fixationends = find(trialevents(:, 2) == scanTargetExitID);
        if ~isempty(fixationstarts) && ~isempty(fixationends)
            [fixationstarts, fixationends] = ...
                lfp_theresa_debounce(fixationstarts, fixationends, trialevents);
        end
        if length(fixationends) == length(fixationstarts) - 1
            fixationends(end+1) = find( ...
                trialevents(:, 2) == scanAllTargetsOffID );
        end
        
        % 'RwdPos':
        taskspecific = [ BDpos(bitand(255, lfp_TrialParams{trial}(rwdPosPN))) ];
        % 'ScanDur':
        scandur = trialevents(find(trialevents(:,2) == scanAllTargetsOffID), 1) ...
            - trialevents(find(trialevents(:,2) == scanAllTargetsOnID), 1);
        taskspecific = [ taskspecific scandur ];
        
        if length(fixationstarts) ~= length(fixationends)
            if length(fixationstarts) > length(fixationends)
                warning('lfp_theresaEventTabulation:fixmismatch1', ...
                    'Skipping trial %d due to unmatched fixation start(s)', ...
                    trial);
            else
                warning('lfp_theresaEventTabulation:fixmismatch2', ...
                    'Skipping trial %d due to unmatched fixation end(s)', ...
                    trial);
            end
            datagood = false;
        elseif length(fixationstarts) == 0
            datagood = false;
        else
            rexnt = trialevents(fixationstarts(1), 1) ...
                - trialevents(find(trialevents(:,2) == scanAllTargetsOnID), 1);
            if isempty(rexnt)
                rexnt = 0;
            end
            
            fixationdurs = trialevents(fixationends, 1) ...
                - trialevents(fixationstarts, 1);
            % 'Pos', 'Dur':
            taskspecific = ...
                [ taskspecific pos(trialevents(fixationstarts(1), 2)) ...
                    fixationdurs(1) ];
            % 'OffDur', 'deltaX', 'deltaY', 'Pos', 'Dur', ...:
            if length(fixationstarts) > 1
                offtargetdurs = trialevents(fixationstarts(2:end), 1) ...
                    - trialevents(fixationends(1:end-1), 1);
                fixationdisps = coords(pos(trialevents( ...
                    fixationstarts(2:end), 2)), :) ...
                    - coords(pos(trialevents( ...
                    fixationstarts(1:end-1), 2)), :) ;
                posns = reshape(...
                    pos(trialevents(fixationstarts(2:end), 2)), ...
                    [], 1);
                taskspecific = ...
                    [ taskspecific reshape( ...
                        [offtargetdurs'
                        fixationdisps(:,1)'
                        fixationdisps(:,2)'
                        posns'
                        fixationdurs(2:end)' ], ...
                        1, [] ) ...
                        ];
                % Auto-indent does not work here, resetting with comment
            end
        end
    end
    if (nargout > 0) && datagood && ( ...
            (fixflag && (trialtype ==  fixationTrial)) ...
            || (scanflag && (trialtype == scanTrial)) )
        table{end+1, 1} = [trial rexnt trialtype trialresult taskspecific];
    elseif (nargout == 0) && datagood && ( ...
            (fixflag && (trialtype ==  fixationTrial)) ...
            || (scanflag && (trialtype == scanTrial)) )
        fprintf(outfid, '%s\t%d\t%d\t%d', ...
            lfp_getTrialID(trial), rexnt, trialtype, trialresult);
        if length(taskspecific) > 252
            fprintf(outfid, '\t%d', taskspecific(1:251));
            fprintf(outfid, '\ttruncated');
        else
            fprintf(outfid, '\t%d', taskspecific);
        end
        fprintf(outfid, '\n');
    end
end
    
if nargout == 0
    fclose(outfid);
end

