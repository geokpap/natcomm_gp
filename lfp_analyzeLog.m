function lfp_analyzeLog(filepath, thresh, varargin)
%lfp_analyzeLog(filepath, thresh)
%lfp_analyzeLog(..., 'LFPtrodeID', LFPtrodeIDstr)
%lfp_analyzeLog(..., 'p2')
%lfp_analyzeLog(..., 'spiketrodeID', spiketrodeIDstr)

%lfp_analyzeLog(filepath, thresh)
% Merged from analyzeLog_bill and lfp_analyzeLog.
% <filepath> is the absolute or relative pathname to the log file (or
% excerpt thereof), thresh is the smallest value of "differential summed
% coherence" that is considered deserving of further analysis.  Checks the
% the last directory in the pathname of the events file to determine
% whether we are starting a new experiment session.
% 1. For each unit, determine if the differential summed coherence
%   exceeds a threshold on at least one of the LFPs; OR, if p2flag,
%   determine whether p2 is below threshold on at least one of the LFPs; if
%   not, create no output.
% 2. For units passing criterion 1, which trode has greatest diff summed
%   coh ("lfp#"), which has smallest p2 ("p2lfp#"), and whether the trode
%   with the max dsc (or min p2 for option 'p2') is the same trode that the
%   unit was on ('TRUE', 'FALSE', or blank if electrode number was not
%   recoverable).
% 3. How many LFPs were above (or below, in case of 'p2') threshold
%
%lfp_analyzeLog(..., 'LFPtrodeID', LFPtrodeIDstr)
% The string of non-digits starting at the beginning of the LFP channel
% name, which is followed by the LFP electrode number in an LFP channel
% name.  Default: 'LFP'.
%lfp_analyzeLog(..., 'p2')
% 'lfp_analyzeLog' style: thresh is interpreted as the significance
% level p2 below which a coherogram is considered deserving of further
% analysis.
%lfp_analyzeLog(..., 'spiketrodeID', spiketrodeIDstr)
% <spiketrodeIDstr> is a string of non-digits the last occurence of which
% is followed by the spike electrode number in a spike channel name.  Any
% value other than 'C' is allowed (e.g. 'SE').  Default is 'T'.

%$Rev: 280 $
%$Date: 2012-07-17 18:42:05 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

global lfp_analyzeLog_output_str;

argnum = 1;
criterion = 'dsc';
LFPtrodeIDstr = 'LFP';
p2flag = false;
spiketrodeIDstr = 'T';
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'LFPtrodeID'
            argnum = argnum + 1;
            LFPtrodeIDstr = varargin{argnum};
        case 'p2'
            p2flag = true;
            criterion = 'p2';
        case 'spiketrodeID'
            argnum = argnum + 1;
            spiketrodeIDstr = varargin{argnum};
        otherwise
            error('lfp_analyzeLog:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

lfp_analyzeLog_output_str = [];
lfp_analyzeLog_cluster_count = 0;
logfid = fopen(filepath);
if logfid == -1
    error('Could not open log file');
end
[pathstr,name] = fileparts(filepath);
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
        sessiondir = fileparts(evtfilepath);
        if ~isequal(prevsessiondir, sessiondir)
            % We have indeed started a new session
            disp(sprintf('Found report for session dir %s at line %d', ...
                sessiondir, linenum));
            if gotdata
                % process last cluster of previous session
                writedata(outfid, prevsessiondir, lfp, spike, dsc, thresh, ...
                    xlim, flim, align, p2, p2flag, spiketrodeIDstr, LFPtrodeIDstr);
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
    elseif length(line) > 46 && ...
            isequal(line(22:46), 'Selected spike channels  ')
        rem = line(47:end);
        nclusts = 0;
        while ~isempty(rem)
            [tok, rem] = strtok(rem);
            if ~isempty(tok)
                nclusts = ...
                    nclusts + 1;
            end
        end
        disp(sprintf('...with analysis of %d clusters', nclusts));
        lfp_analyzeLog_cluster_count = ...
            lfp_analyzeLog_cluster_count + nclusts;
    elseif length(line) > 31 && isequal(line(22:32), 'Coherogram:')
        [lfpname, rem] = strtok(line(33:end));
        if ~( length(lfpname) >= length(LFPtrodeIDstr) ...
                && isequal(lfpname(1:length(LFPtrodeIDstr)), LFPtrodeIDstr) )
            warning('lfp_analyzeLog:e4', ...
                'Bad LFP name "%s"', lfpname );
            continue;
        end
        lfp{end+1} = lfpname;
        [junk, rem] = strtok(rem);
        [spike{end+1}, rem] = strtok(rem);
        if (length(spike) > 1) && ~isequal(spike{end}, spike{end-1})
            % We just started a new cluster, process the one that's
            % finished
            writedata(outfid, sessiondir, lfp(1:end-1), ...
                spike(1:end-1), dsc, thresh, ...
                xlim, flim, align, p2, p2flag, spiketrodeIDstr, LFPtrodeIDstr);
            % And remove the previous cluster from the data lists
            lfp = lfp(end);
            spike = spike(end);
            dsc = [];
            p2 = [];
        end
        if ~isempty(rem)
            warning('lfp_analyzeLog:e1', ...
                'Trouble parsing "%s"', line );
        end
    elseif length(line) > 53 && isequal(line(22:50), 'differential summed coherence')
        dsc(end+1) = sscanf(line(54:end), '%f');
        gotdata = true;
    elseif length(line) > 67 && isequal(line(22:64), 'lfp_spec: significance level of cell counts')
        p2(end+1) = sscanf(line(68:end), '%f');
        gotdata = true;
    elseif length(line) > 42 && isequal(line(22:43), 'Error while processing')
        warning('lfp_analyzeLog:e2', ...
            'Logged error: "%s"', line );
    end
end
% We have reached EOF
if gotdata
    writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
        xlim, flim, align, p2, p2flag, spiketrodeIDstr, LFPtrodeIDstr);
end
fclose(logfid);
lfp_analyzeLog_output_str = [
    sprintf('logfile:\t%s\tthresh:\t%d\tclusters:\t%d\tcriterion:\t%s\n\n', ...
    filepath, thresh, lfp_analyzeLog_cluster_count, criterion  ) ...
    sprintf('sessiondir\tcluster\tmaxdsc\tlfp#\tsametrode\t# trodes\tXLim\tFreqLim\tAlign\tp2lfp#\tminp2\n') ...
    lfp_analyzeLog_output_str ];
fprintf(outfid, '%s', lfp_analyzeLog_output_str);
fclose(outfid);


function writedata(outfid, sessiondir, lfp, spike, dsc, thresh, ...
    xlim, flim, align, p2, p2flag, spiketrodeIDstr, LFPtrodeIDstr)
global lfp_analyzeLog_output_str;
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
        shortones = [shortones ' p2'];
    end
    warning('lfp_analyzeLog:e3', ...
        'Incomplete data: %s', shortones );
end
suprathresh = find(dsc>=thresh);
subthresh = find(p2<thresh);
if p2flag
    numtrodes = length(subthresh);
else
    numtrodes = length(suprathresh);
end
if (p2flag && ~isempty(subthresh)) || (~p2flag && ~isempty(suprathresh))
    [maxdsc, idx] = max(dsc);
    [minp2, idxp2] = min(p2);
    if ~isempty(idx)
        lfptrodenum = str2num(lfp{idx}(length(LFPtrodeIDstr)+1:end));
    else
        lfptrodenum = [];
    end
    if ~isempty(idxp2)
        p2lfptrodenum = str2num(lfp{idxp2}(length(LFPtrodeIDstr)+1:end));
    else
        p2lfptrodenum = [];
    end
    % In spike channel names, cluster numbers are always indicated by 'C'
    % followed by digits followed by end-of-string.  This means that the
    % cluster number always follows the last 'C' in the string.  However,
    % spike electrode numbers are heterogeneously identified by an
    % arbitrary non-digit string <spiketrodeIDstr> that may occur in the
    % middle of the channel name (known values at this point are 'T' and
    % 'SE'; anything other than 'C' will be fine).
    spiketrodeIDidx = strfind(spike{idx}, spiketrodeIDstr);
    Cidx = find(spike{idx} == 'C');
    sametrode = '';
    if ~isempty(spiketrodeIDidx) && ~isempty(Cidx)
        spiketrodenum = str2num(spike{idx}(...
            spiketrodeIDidx(end)+length(spiketrodeIDstr) : Cidx(end)-1 ));
        if (~p2flag && ~isempty(lfptrodenum))
            sametrode = (lfptrodenum == spiketrodenum);
        elseif (p2flag && ~isempty(p2lfptrodenum))
            sametrode = (p2lfptrodenum == spiketrodenum);
        end
    else
        warning('lfp_analyzeLog:badname', ...
            'Could not parse spike channel name "%s"', spike{idx} );
    end
    % The terminating newline must go in a separate sprintf to handle the
    % case where data are missing
    lfp_analyzeLog_output_str = [ lfp_analyzeLog_output_str ...
        sprintf('%s\t%s\t%f\t%i\t%s\t%d\t%s\t%s\t%s\t%d\t%d', ...
        sessiondir, spike{idx}, maxdsc, lfptrodenum, mat2str(sametrode), ...
        numtrodes, mat2str(xlim'), mat2str(flim'), ...
        dg_thing2str(align), p2lfptrodenum, minp2 ) ...
        sprintf('\n') ];
end
