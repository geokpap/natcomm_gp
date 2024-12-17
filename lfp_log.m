function lfp_log(message)
%lfp_log(message)
% If lfp_LogMode is empty, then each 

%$Rev: 256 $
%$Date: 2012-01-30 17:38:09 -0500 (Mon, 30 Jan 2012) $
%$Author: dgibson $

global lfp_LogFileName

if isequal(lfp_LogFileName, 'none')
    return
end
if isempty(lfp_LogFileName)
    error('lfp_LogFileName is empty.');
end
files = dir(lfp_LogFileName);

if ~isempty(files)
    % Truncate head of file if necessary
    if files.bytes > 10e6
        fid = fopen(lfp_LogFileName, 'r+');
        if fid == -1
            error('lfp_log:open1', ...
                'Could not open log file %s', lfp_LogFileName );
            return
        end
        % skip first half of file
        status = fseek(fid, fix(files.bytes/2), 'bof');
        if status == -1
            warning('lfp_log:badseek', ...
                'Seek failed in file %s', lfp_LogFileName );
            fclose(fid);
        else
            % read rest of file, find the end of the first (usually incomplete)
            % line of text
            chars = fread(fid, inf, 'uchar');
            fclose(fid);
            fid = fopen(lfp_LogFileName, 'w');
            if fid == -1
                error('lfp_log:open3', ...
                    'Could not open log file %s for writing', lfp_LogFileName );
                return
            end
            % The PC eol ('\r\n') and the Unix eol ('\n') both end in '\n',
            % so regardless of what type of eols we have, we want to
            % truncate the head through the first '\n'.  Note that the two
            % different styles might be mixed together on a PC because of
            % the large number of sprintf('blabla\n') formulations in
            % lfp_lib.
            eol = sprintf('\n');
            eolidx = 1;
            while eolidx <= length(chars) && ~isequal( ...
                    chars(eolidx - length(eol) + 1 : eolidx), eol )
                eolidx = eolidx + 1;
            end
            if  length(chars) >= eolidx
                count = fwrite(fid, chars(eolidx+1:end), 'uchar');
                if count ~= length(chars) - eolidx
                    error('lfp_log:trunc', ...
                        'Failed while writing truncated log file.');
                end
            end
            fclose(fid);
        end
    end
end

fid = fopen(lfp_LogFileName, 'a');
for numattempts = 1:1000
    if fid == -1
        if mod(numattempts, 10) == 0
            warning('lfp_log:open2', ...
                'Could not open log file %s', lfp_LogFileName );
        end
        fid = fopen(lfp_LogFileName, 'a');
    else
        break
    end
end
if fid == -1
    error('lfp_log:open4', ...
        'Could not open log file %s with %d attempts', ...
        lfp_LogFileName, numattempts );
    return
end
fprintf(fid, '%s %s\n', datestr(now, 0), message);
status = fclose(fid);
if status == -1
    error('lfp_log:close', ...
        'Could not close log file fid %d, %s', fid, lfp_LogFileName );
    return
end
