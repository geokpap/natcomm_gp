function analyzeLog_bill(filepath, thresh)
% For analyzing log output from bulkProcess_bill(4, ...).  <filepath> is the
% absolute or relative pathname to the log file (or excerpt thereof),
% thresh is the smallest value of "differential summed coherence" that is
% considered deserving of further analysis.  See function writedata below
% for details.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

logfid = fopen(filepath);
if logfid == -1
    error('Could not open log file');
end
[pathstr,name,ext,versn] = fileparts(filepath);
OutputFileName = [name '.xls'];
[OutputFileName, OutputPathName] = uiputfile( ...
    fullfile(pathstr, OutputFileName) );
if isequal(OutputFileName, 0)
    fclose(logfid);
    return
end
outfid = fopen(fullfile(OutputPathName, OutputFileName), 'w');
if outfid == -1
    fclose(logfid);
    error('Could not open output file');
end
fprintf(outfid, 'sessiondir\tcluster\tmaxdsc\tlfp#\tsametrode\t# trodes\tXLim\tFreqLim\tAlign\n');

sessiondir = '';
line = '';
linenum = 0;
gotdata = false;
xlim = [];
flim = [];
align = [];
lfp = {};
spike = {};
dsc = [];
while ~isequal(line, -1)
    line = fgetl(logfid);
    linenum = linenum + 1;
    if gotdata && length(line) > 2 ...
            && isequal(line(3:end), 'STARTING FROM CLEARED MEMORY')
        % We just started another session, process last cluster of last
        % session
        writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
            xlim, flim, align);
        gotdata = false;
    elseif length(line) > 36 && isequal(line(22:37), 'Read events file')
        evtfilepath = line(40:end-1);
        [sessiondir,name,ext,versn] = fileparts(evtfilepath);
        disp(sprintf('Found report for session dir %s at line %d', ...
            sessiondir, linenum));
    elseif length(line) > 12 && isequal(line(2:13), 'lfp_XLimAll=')
        stuff = sscanf(line(2:end), ...
            'lfp_XLimAll=[%f %f], lfp_FreqLim=[%f %f], lfp_AlignmentRef=%i, lfp_CLimAll=[%f %f]' );
        if length(stuff) >1
            xlim = stuff(1:2);
        else
            xlim = [];
        end
        if length(stuff) > 3
            flim = stuff(3:4);
        else
            flim = [];
        end
        alignidx = strfind(line, 'lfp_AlignmentRef=');
        if isempty(alignidx)
            align = [];    
        else
            alignidx = alignidx(1);
            commaidx = strfind(line(alignidx:end), ',');
            if ~isempty(commaidx)
                commaidx = commaidx(1) + alignidx - 1;
            else
                commaidx = length(line) + 1;
            end
            alignstr = line(alignidx + 17 : commaidx - 1);
            if alignstr(1) == '['
                alignstr = alignstr(2:end-1);
            end
            align = reshape(sscanf(alignstr, '%i'), 1, []);
        end
        lfp = {};
        spike = {};
        dsc = [];
    elseif length(line) > 31 && isequal(line(22:32), 'Coherogram:')
        [lfp{end+1}, rem] = strtok(line(33:end));
        [junk, rem] = strtok(rem);
        [spike{end+1}, rem] = strtok(rem);
        if (length(spike) > 1) && ~isequal(spike{end}, spike{end-1})
            % We just started a new cluster, process the one that's
            % finished
            writedata(outfid, sessiondir, lfp(1:end-1), ...
                spike(1:end-1), dsc, thresh, ...
                xlim, flim, align);
            % And remove the previous cluster from the data lists
            lfp = lfp(end);
            spike = spike(end);
            dsc = [];
        end
        if ~isempty(rem)
            warning('analyzeLog_bill:e1', ...
                'Trouble parsing "%s"', line );
        end
    elseif length(line) > 49 && isequal(line(22:50), 'differential summed coherence')
        dsc(end+1) = sscanf(line(54:end), '%f');
        gotdata = true;
    elseif length(line) > 42 && isequal(line(22:43), 'Error while processing')
        warning('analyzeLog_bill:e2', ...
            'Logged error: "%s"', line );
    end
end
if gotdata
    writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
        xlim, flim, align);
end
fclose(logfid);
fclose(outfid);


function writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
    xlim, flim, align)
%1. Make a list of units where the differential summed coherence
%   exceeds a threshold on at least one of the LFPs 
%2. For units passing criterion 1, which trode has greatest diff summed
%   coh, and whether that is the same trode that the unit was on 
%3. How many LFPs were above threshold
%
% We assume that every LFP file has a name consisting of 3 non-digits
% followed a string of digits denoting the electrode number.
if (length(lfp) ~= length(spike)) || (length(lfp) ~= length(dsc))
    maxlen = max([length(lfp) length(spike) length(dsc)]);
    shortones = '';
    if length(lfp) < maxlen
        lfp(end+1:maxlen) = {''};
        shortones = [shortones ' lfp'];
    end
    if length(spike) < maxlen
        spike(end+1:maxlen) = {''};
        shortones = [shortones ' spike'];
    end
    if length(dsc) < maxlen
        dsc(end+1:maxlen) = -1;
        shortones = [shortones ' dsc'];
    end
    warning('analyzeLog_bill:e3', ...
        'Incomplete data: %s', shortones );
end
suprathresh = find(dsc>=thresh);
if ~isempty(suprathresh)
    [maxdsc, idx] = max(dsc);
    lfptrodenum = str2num(lfp{idx}(4:end));
    Tidx = find(spike{idx} == 'T');
    Cidx = find(spike{idx} == 'C');
    if ~isempty(Tidx) && ~isempty(Cidx)
        spiketrodenum = str2num(spike{idx}(Tidx(end)+1:Cidx(end)-1));
        sametrode = (lfptrodenum == spiketrodenum);
    else
        sametrode = '';
    end
    fprintf(outfid, '%s\t%s\t%f\t%i\t%s\t%d\t%s\t%s\t%s', ...
        sessiondir, spike{idx}, maxdsc, lfptrodenum, mat2str(sametrode), ...
        length(suprathresh), mat2str(xlim'), mat2str(flim'), dg_thing2str(align) );
    fprintf(outfid, '\n');
end
