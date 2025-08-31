function numfiles = dg_filecat(infiles, outfile, verify)
%dg_filecat(infiles, outfile) concatenates the contents of one or more
% input files to an output file.

% If an input file cannot be read, a warning is issued and the file is
% skipped.  If the output file already exists or cannot be written, an
% error is raised.
% INPUTS
%   infiles: a cell array of strings, each of which points to a file to be
%       read and copied to the output.
%   outfile: a string that points to the destination file to which output
%       is copied.
%   verify: may be 0, 1, or 2.
%       0:  no verification is done and the files are simply concatenated
%       including any header lines they might contain. 
%       1:  the first line of the first file is compared to the first line
%       of each subsequent file, and if they are not identical, then the
%       file is skipped and a warning issued.  (This can be used to prevent
%       tab-delimited text spreadsheet files from being concatenated when
%       the headers don't match.)  Also, the first line is only copied from
%       the first file; subsequent first lines are omitted.
%       2:  similar to 1, but instead of using the first file's header line
%       as the reference, all the files are scanned and the longest header
%       line is used as the reference; if there are any headers that
%       containing column names that mismatch the ones in the reference
%       header, an error is raised.  
% NOTES
%   This program is written to work with both DOS-style ('\r\n') and
%   Unix-style ('\n') line terminators; however, it will not work with
%   <verify> set if different files use different styles of terminators for
%   the first line.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 3
    error('Call as dg_filecat(infiles, outfile, verify)');
end
numfiles = 0;
if exist(outfile)
    error('The output file already exists.');
end
outfid = fopen(outfile, 'w');
if outfid == -1
    error('Could not open output file for writing.');
end
switch verify
    case 0
        startbyte = 1;
    case 1
        infid = fopen(infiles{1}, 'r');
        if infid == -1
            fclose(outfid);
            error('Could not open first input file.');
        end
        refheader = reshape(fgets(infid), 1, []);
        fclose(infid);
        startbyte = length(refheader) + 1;
        fprintf(outfid, '%s', refheader);
    case 2
        for k = 1:numel(infiles)
            infid = fopen(infiles{k}, 'r');
            if infid == -1
                fclose(outfid);
                error('Could not open input file "%s".', infiles{k});
            end
            headers{k} = reshape(fgets(infid), 1, []);
            headerlengths(k) = length(headers{k});
            fclose(infid);
        end
        [C, I] = max(headerlengths);
        refheader = headers{I};
        if length(refheader) > 1 && refheader(end-1) == 13
            EOLlen = 2;
        else
            EOLlen = 1;
        end
        for k = 1:numel(infiles)
            if ~strncmp(headers{k}, refheader, length(headers{k}) - EOLlen)
                fclose(outfid);
                error('Header mismatch in file "%s".', infiles{k});
            end
        end
        fprintf(outfid, '%s', refheader);
end
for k = 1:numel(infiles)
    infid = fopen(infiles{k}, 'r');
    if infid == -1
        warning('Could not open input file "%s".', infiles{k});
        continue
    end
    slurp = fread(infid, Inf, 'uchar');
    fclose(infid);
    if (verify == 1) && ~isequal(reshape(slurp(1:length(refheader)), 1, []), refheader)
        warning('Header mismatch, skipping file "%s"', infiles{k});
        continue
    end
    numfiles = numfiles + 1;
    if verify == 2
        startbyte = length(headers{k}) + 1;
    end
    fwrite(outfid, slurp(startbyte:end), 'uchar');
end
fclose(outfid);
