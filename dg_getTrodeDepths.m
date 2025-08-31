function [depths, sessionIDs] = dg_getTrodeDepths( ...
    dataroots, trodenum, dirformat)
%[depths, sessionIDs] = dg_getTrodeDepths(dataroots, trodenum, dirformat)
% Crawls all directory trees listed in <dataroots> in search of
% session-by-session depths for electrode <trodenum>.  Works for
% <dirformat> = 'ken'.
%INPUTS
% dataroots: cell string array
% trodenum: integer scalar
% dirformat: string
%OUTPUTS
% depths: numeric array with one row per session and at least 2 columns. 
%   Col 1 = datenum for session without hrs, min or sec (NaN if session ID
%   is not convertible), col 2 = depth (units and reference point depend on
%   <dirformat>).  May contain -1 to mark a session that exists but does
%   not contain valid depth data.  Additional columns for:
%       ken data: 
%           col 3 = stim CT
%           col 4 = stim AN
% sessionIDs: cell string vector with one row per session.

%$Rev: 169 $
%$Date: 2013-03-01 17:48:02 -0500 (Fri, 01 Mar 2013) $
%$Author: dgibson $

depths = zeros(0, 2);
sessionIDs = {};

if ~isequal(dirformat, 'ken')
    error('dg_getTrodeDepths:dirformat', ...
        'Not available for dirformat "%s".', dirformat);
end

for diridx = 1:length(dataroots)
    files = dir(dataroots{diridx});
    is_dir = cell2mat({files.isdir});
    mysessionIDs = {files(is_dir).name};
    for sessidx = 1:length(mysessionIDs)
        dirstr = mysessionIDs{sessidx};
        toks = regexp(dirstr, '^(\d+)-(\d+)-(\d+).*$', 'tokens');
        if isempty(toks)
            datevector = NaN;
        else
            datevector = [str2double(toks{1}) 0 0 0];
        end
        if ismember(dirstr, {'.' '..'})
            continue
        end
        % Extract depth data
        if exist(fullfile(dataroots{diridx}, dirstr, 'aloc1.txt'), 'file')
            thecmd = sprintf( 'grep UN:%d ''%s/aloc1.txt''', ...
                trodenum, fullfile(dataroots{diridx}, dirstr) );
            [s,r]=system(thecmd);
            switch s
                case 0
                    % fine, do nothing
                case 1
                    % no line for that electrode
                    depthstr = '-1';
                otherwise
                    error('oops')
            end
            depthtok = regexpi(r(1:end-2), 'DP:([\d.]+)$', 'tokens');
            if isempty(depthtok) || isempty(depthtok{1})
                depthstr = '-1';
            else
                depthstr = depthtok{1}{1};
            end
        else
            depthstr = '-1';
        end
        % Extract stim data
        if exist(fullfile(dataroots{diridx}, dirstr, 'stim1.txt'), 'file')
            thecmd = sprintf( 'grep CT:%d ''%s/aloc1.txt''', ...
                trodenum, fullfile(dataroots{diridx}, dirstr) );
            [s,r]=system(thecmd);
            switch s
                case 0
                    % fine, do nothing
                case 1
                    % no line for that electrode
                    CTstr = '-1';
                otherwise
                    error('barf')
            end
            stimtok = regexpi(r(1:end-2), 'CT:([\d.]+),.*$', 'tokens');
            if isempty(stimtok) || isempty(stimtok{1})
                CTstr = '-1';
            else
                CTstr = stimtok{1}{1};
            end
        else
            CTstr = '-1';
        end
        if exist(fullfile(dataroots{diridx}, dirstr, 'stim1.txt'), 'file')
            thecmd = sprintf( 'grep AN:%d ''%s/aloc1.txt''', ...
                trodenum, fullfile(dataroots{diridx}, dirstr) );
            [s,r]=system(thecmd);
            switch s
                case 0
                    % fine, do nothing
                case 1
                    % no line for that electrode
                    ANstr = '-1';
                otherwise
                    error('barf')
            end
            stimtok = regexpi(r(1:end-2), 'CT:([\d.]+),.*$', 'tokens');
            if isempty(stimtok) || isempty(stimtok{1})
                ANstr = '-1';
            else
                ANstr = stimtok{1}{1};
            end
        else
            ANstr = '-1';
        end
        try
            depths(end+1, 1) = datenum(datevector); %#ok<*AGROW>
        catch
            depths(end+1, 1) = NaN;
        end
        depths(end, 2) = str2double(depthstr);
        depths(end, 3) = str2double(CTstr);
        depths(end, 4) = str2double(ANstr);
        sessionIDs{end+1,1} = dirstr;
    end
end


