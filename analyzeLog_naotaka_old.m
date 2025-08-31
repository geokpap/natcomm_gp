function analyzeLog_naotaka(filepath, thresh)
%analyzeLog_naotaka(filepath, thresh)
% Based on analyzeLog_bill.  <filepath> is the
% absolute or relative pathname to the log file (or excerpt thereof),
% thresh is the significance level p2 below which a coherogram is
% considered deserving of further analysis.  See function writedata below
% for details.
%
% Modified 3/2/2005 by DG to report all clusters with signif level of
% significant cell counts p2<thresh, and two additional columns containing
% LFP channel with smallest p2 and that p2 value.  Since these analyses
% require multiple lfp_read's to process one experiment session, this
% script also checks the the last directory in the pathname of the events
% file to determine whether we are starting a new experiment session.

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
fprintf(outfid, 'sessiondir\tcluster\tmaxdsc\tlfp#\tsametrode\t# trodes\tXLim\tFreqLim\tAlign\tp2lfp#\tminp2\n');

sessiondir = '';
line = '';
linenum = 0;
gotdata = false;
readflag = false;
xlim = [];
flim = [];
align = [];
lfp = {};
spike = {};
dsc = [];
while ~isequal(line, -1)
    line = fgetl(logfid);
    linenum = linenum + 1;
    if length(line) > 2 ...
            && isequal(line(3:end), 'STARTING FROM CLEARED MEMORY')
        % We may have just started another session, remember to check
        readflag = true;
    elseif readflag && length(line) > 36 ...
            && isequal(line(22:37), 'Read events file')
        readflag = false;
        %determine if this is actually a new session
        evtfilepath = line(40:end-1);
        prevsessiondir = sessiondir;
        [sessiondir,name,ext,versn] = fileparts(evtfilepath);
        if ~isequal(prevsessiondir, sessiondir)
            % We have indeed started a new session
            disp(sprintf('Found report for session dir %s at line %d', ...
                sessiondir, linenum));
            if gotdata
                % process last cluster of previous session
                writedata(outfid, prevsessiondir, lfp, spike, dsc, thresh, ...
                    xlim, flim, align, p2);
                gotdata = false;
            end
        end
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
        p2 = [];
    elseif length(line) > 31 && isequal(line(22:32), 'Coherogram:')
        [cscname, rem] = strtok(line(33:end));
        if ~isequal(cscname(1:3), 'CSC')
            warning('analyzeLog_naotaka:e4', ...
                'Bad CSC name "%s"', cscname(1:3) );
            continue;
        end
        lfp{end+1} = cscname;
        [junk, rem] = strtok(rem);
        [spike{end+1}, rem] = strtok(rem);
        if (length(spike) > 1) && ~isequal(spike{end}, spike{end-1})
            % We just started a new cluster, process the one that's
            % finished
            writedata(outfid, sessiondir, lfp(1:end-1), ...
                spike(1:end-1), dsc, thresh, ...
                xlim, flim, align, p2);
            % And remove the previous cluster from the data lists
            lfp = lfp(end);
            spike = spike(end);
            dsc = [];
            p2 = [];
        end
        if ~isempty(rem)
            warning('analyzeLog_naotaka:e1', ...
                'Trouble parsing "%s"', line );
        end
    elseif length(line) > 53 && isequal(line(22:50), 'differential summed coherence')
        dsc(end+1) = sscanf(line(54:end), '%f');
        gotdata = true;
    elseif length(line) > 67 && isequal(line(22:64), 'lfp_spec: significance level of cell counts')
        p2(end+1) = sscanf(line(68:end), '%f');
        gotdata = true;
    elseif length(line) > 42 && isequal(line(22:43), 'Error while processing')
        warning('analyzeLog_naotaka:e2', ...
            'Logged error: "%s"', line );
    end
end
% We have reached EOF
if gotdata
    writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
        xlim, flim, align, p2);
end
fclose(logfid);
fclose(outfid);


function writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
    xlim, flim, align, p2)
%1. Determine whether p2 is below threshold on at least
%   one of the LFPs; if not, print no output.
%2. For units passing criterion 1, which trode has greatest diff summed
%   coh, and whether that is the same trode that the unit was on (1 = yes,
%   blank = no)
%3. How many LFPs were above threshold
%
% We assume that every LFP file has a name consisting of 3 non-digits
% followed a string of digits denoting the electrode number.
if (length(lfp) ~= length(spike)) ...
        || (length(lfp) ~= length(dsc)) ...
        || (length(lfp) ~= length(p2))
    maxlen = max([length(lfp) length(spike) length(dsc) length(p2)]);
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
    if length(p2) < maxlen
        dsc(end+1:maxlen) = -1;
        shortones = [shortones ' dsc'];
    end
    warning('analyzeLog_naotaka:e3', ...
        'Incomplete data: %s', shortones );
end
subthresh = find(p2<thresh);
if ~isempty(subthresh)
    [maxdsc, idx] = max(dsc);
    [minp2, idxp2] = min(p2);
    lfptrodenum = str2num(lfp{idx}(4:end));
    p2lfptrodenum = str2num(lfp{idxp2}(4:end));
    Eidx = find(spike{idx} == 'E');
    Cidx = find(spike{idx} == 'C');
    if ~isempty(Eidx) && ~isempty(Cidx)
        spiketrodenum = str2num(spike{idx}(Eidx(end)+1:Cidx(end)-1));
        sametrode = (lfptrodenum == spiketrodenum);
    else
        sametrode = '';
    end
    fprintf(outfid, '%s\t%s\t%f\t%i\t%s\t%d\t%s\t%s\t%s\t%d\t%d', ...
        sessiondir, spike{idx}, maxdsc, lfptrodenum, mat2str(sametrode), ...
        length(subthresh), mat2str(xlim'), mat2str(flim'), ...
        dg_thing2str(align), p2lfptrodenum, minp2 );
    fprintf(outfid, '\n');
end
