function dg_AM2Mat(pathname, varargin)
% Convert Ayano format multichannel photometry data with Arduino task
% control (*.csv files) to lfp_lib-compatible Matlab format.

%INPUTS
% pathname: absolute or relative pathname to one file to convert.
%OUTPUTS
% All outputs are sent to output files.  An Arduino events file produces
% one '.mat' CSC format output file with the same base filename as the
% input; a photometry file produces one lfp_lib '.evtsav' for each channel,
% named by concatenating the column header onto the base filename of
% <pathname>.
%OPTIONS
% 'dest', destdir - writes output to <destdir> instead of to the same
%   directory as the source file.
%NOTES
% Edited from a copy of dg_Nlx2Mat.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

destdir = '';
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dest'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        otherwise
            error('dg_Nlx2Mat:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[pathstr,name,ext] = fileparts(pathname);
if ~isequal(lower(ext), '.csv')
    error('dg_AM2Mat:input', ...
        'Input file must have ".csv" file extension.');
end
if isempty(destdir)
    destdir = pathstr;
end

% The only way to tell what kind of file it is is by reading it.  The two
% possibilities are multichannel photometry data or Arduino task events.
% Attempt to read CSV file:
delim = sprintf(',');
newline = sprintf('\n'); %#ok<SPRINTFN>
bytes = fileread(pathname);
lines = strsplit(bytes, newline, 'CollapseDelimiters', false);
cells = cell(length(lines), 0);
for linenum = 1:length(lines)
    values = strsplit(lines{linenum}, delim, 'CollapseDelimiters', false);
    cells(linenum, 1:length(values)) = values;
end
numrows = size(cells, 1);
numcols = size(cells, 2);
if numrows == 0 || numcols < 2
    error('dg_AM2Mat:input', ...
        'Input file must be non-empty and in CSV format.');
end

switch numcols
    case 3
        fprintf('Read Arduino events file.\n');
        % Arduino events file; reformat and save as lfp_lib '.evtsav'.
        % Column headers are optional:
        if all(isnan(str2double(cells(1,:))))
            cells(1,:) = [];
        end
        % Trailing empty rows are prohibited:
        while all(isnan(str2double(cells(end,:))))
            cells(end,:) = [];
        end
        ishandled = false(size(cells, 1), 1);
        allTS = str2double(cells(:,1)); % in ms
        allTTL = str2double(cells(:, 2));
        allval = str2double(cells(:, 3));
        evtTTLs = [ 7 8 13 14 16 46 47 ];
        isevent = ismember(allTTL, evtTTLs);
        TS = allTS(isevent) * 1e-3; % convert ms to s
        TTL = allTTL(isevent);
        lfp_save_events = [TS TTL];
        ishandled(isevent) = true;
        % In addition to classic lfp_TrialParams, Ayano also uses
        % session-wide parameters that are worth saving, although lfp_lib
        % can't actually use them.  There should be exactly one of each of
        % these in the entire file.  They are indexed by TTL ID, so most of
        % them will be empty, and each entry is the 'behavioralDataValue'
        % for the indexed ID.
        sessTTLs = [1:4 11 12 40 49 56 59 68 92 93 116 120:134];
        isSessParam = ismember(allTTL, sessTTLs);
        if length(unique(allTTL(isSessParam))) ~= sum(isSessParam)
            warning('dg_AM2Mat:sessTTLs', ...
                'Extra sessTTLs.');
        end
        % Session parameters must have the exact same timestamp as the TTL
        % ID=1.
        evt1TS = allTS(allTTL == 1);
        sessParamidx = find(allTS == evt1TS);
        if length(sessParamidx) < sum(isSessParam)
            warning('dg_AM2Mat:sessParams', ...
                'There are displaced session parameters.');
        end
        AM_sess_params(allTTL(sessParamidx)) = allval(sessParamidx);
        ishandled(sessParamidx) = true;
        % All other events are presumed to be trial parameters.  Trial
        % parameters must have the exact same timestamp as the event 7
        % ('TTL_EVT.TRIAL_START') that identifies the trial.
        evt7idx = find(allTTL == 7);
        numtrials = length(evt7idx);
        lfp_save_params = cell(numtrials, 1);
        for trialidx = 1:numtrials
            trialTS = allTS(evt7idx(trialidx));
            TPidx = find(allTS == trialTS);
            if length(TPidx) ~= 7
                warning('dg_AM2Mat:TPidx', ...
                    'Wrong number of params for trial index %d.', ...
                    trialidx);
            else
                lfp_save_params{trialidx, 1} = ...
                    reshape(allval(TPidx), 1, []); 
                ishandled(TPidx) = true;
            end
        end
        if any(~ishandled)
            warning('dg_AM2Mat:unhandled', ...
                'There are unhandled event records.');
        end
        save(fullfile(destdir, [name '.evtsav']), ...
            'lfp_save_events', 'lfp_save_params', 'AM_sess_params');
    otherwise
        % Multichannel photometry, may contain any number of channels.
        % Each channel goes in its own separate '.mat' format CSC file.
        numch = numcols - 1;
        fprintf('Read photometry file with %d channels.\n', numch);
        % Column headers are required:
        if ~all(isnan(str2double(cells(1,:))))
            error('dg_AM2Mat:hdr', ...
                'Column headers are required in photometry files.');
        end
        hdr = cells(1, 2:end);
        for chnum = 1:numch
            dg_Nlx2Mat_Timestamps = reshape( ...
                str2double(cells(2:end, 1)), 1, [] );
            dg_Nlx2Mat_Samples = reshape( ...
                str2double(cells(2:end, chnum + 1)), 1, [] );
            dg_Nlx2Mat_SamplesUnits = 'arbs';
            matfilename = fullfile(destdir, [ name hdr{chnum} '.mat' ]);
            save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
                'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
        end
end