function lfp_makeCohFigs(filepath, varargin) 
%lfp_makeCohFigs(filepath)
% Reads analyzeLog_bill output file in tab-delimited text format from
% <filepath>, and for each row of the file, creates two coherograms:  one
% with [nw,k] = [2,3] and one with [nw,k] = [5,9] (both with minspikes =
% 2*k).  Uses the values in the "sessiondir", "cluster", "lfp#" and "Align"
% columns to locate the data according to rat directory structure
% conventions, and set lfp_AlignmentRef.  Reads spike data from *.Tnn files
% (will not work with *.TTn files).  Uses previously set values of
% lfp_XLimAll, lfp_FreqLim, lfp_CLimAll, setup name, etc. Runs each
% coherogram in its own try...catch block for error handling. Saves each
% figure to a Matlab *.fig file in the current Matlab working directory
% with name <clustername><lfpname>_<events>_<nw>-<k> where <clustername> is
% the cluster's name in lfp_SpikeNames, <lfpname> is the LFP file name in
% lfp_FileNames, <events> is the value of lfp_AlignmentRef with multiple
% values separated by '-'. 
%lfp_makeCohFigs(..., 'LFPpre', prefix)
%lfp_makeCohFigs(..., 'LFPext', ext)
% Constructs LFP file names as <prefix><lfpnum>.<ext> where <lfpnum>
% is read from the input file and <ext> is defined in the lfp_setup
% file.  Defaults: prefix='LFP', ext='DAT'.
%lfp_makeCohFigs(..., 'monkey')
% Locates data according to monkey directory structure conventions and
% reads spike data from Naotaka cluster files.
%lfp_makeCohFigs(..., 'nwk', nwk)
% Creates one coherogram for each row of nwk.  First column represents
% values of nw, second column is values of k.
%lfp_makeCohFigs(..., 'p2')
% Reads analyzeLog_naotaka output and uses the "p2lfp#" column instead of
% the "lfp#" column to choose the LFP file.
%lfp_makeCohFigs(..., 'pad', padfactor)
% Sets lfp_spec pad factor; default 2.
%lfp_makeCohFigs(..., 'trialselect', evalstr)
% After loading each spike file, <evalstr> is evaluated, the intention
% being to provide a flexible way of selecting trials.  However, <evalstr>
% may be any Matlab-executable string.  If this option is not given, or if
% <evalstr> is empty, the default action is
%         lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
%         lfp_selectByDuration('and');
%lfp_makeCohFigs(..., 'win', moving_win)
% Sets value of <moving_win> used in call to lfp_spec; default value in
% this function is [1 0.1].

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

argnum = 1;
evalstr = '';
LFPpre = 'LFP';
LFPext = 'DAT';
monkeyflag = false;
moving_win = [1 0.1];
p2flag = false;
padfactor = 2;
nwk = [2 3; 5 9];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'LFPext'
            argnum = argnum + 1;
            LFPext = varargin{argnum};
        case 'LFPpre'
            argnum = argnum + 1;
            LFPpre = varargin{argnum};
        case 'monkey'
            monkeyflag = true;
        case 'nwk'
            argnum = argnum + 1;
            nwk = varargin{argnum};
        case 'p2'
            p2flag = true;
        case 'pad'
            argnum = argnum + 1;
            padfactor =  varargin{argnum};
        case 'trialselect'
            argnum = argnum + 1;
            evalstr = varargin{argnum};
        case 'win'
            argnum = argnum + 1;
            moving_win = varargin{argnum};
        otherwise
            error('lfp_makeCohFigs:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~(isequal(size(moving_win), [1 2])) ...
        && isequal(class(moving_win), 'double')
    error('lfp_makeCohFigs:badwin', ...
        '<moving_win> must be a 1x2 array of type double.' );
end

if ~( size(nwk, 2) == 2 ...
        && isequal(class(nwk), 'double') )
    error('lfp_makeCohFigs:badnwk1', ...
        '<nwk> must be a 2-column array of type double.' );
end
if ~isequal(round(nwk(:,2)), nwk(:,2))
    error('lfp_makeCohFigs:badnwk2', ...
        '<k> values must be integers.' );
end
if any(nwk(:,2) > 2*nwk(:,1))
    warning('lfp_makeCohFigs:badnwk3', ...
        '<k> values should be less than 2*<nw>.' );
end

fid = fopen(filepath);
if fid == -1
    error('Could not open log file');
end

%
% Read the file
%
line = '';
linenum = 0;
if p2flag
    lfpcol = 10;
else
    lfpcol = 4;
end
sessioncol = 1;
clustercol = 2;
aligncol = 9;
data = [];
while 1
    line = fgetl(fid);
    % find EOF:
    if isequal(line, -1)
        break
    end
    % skip header lines:
    if isempty(line) ...
            || (length(line) > 6 && isequal(line(1:7), 'logfile')) ...
            || (length(line) > 9 && isequal(line(1:10), 'sessiondir'))
        continue
    end
    % collect data:
    linenum = linenum + 1;
    delims = strfind(line, sprintf('\t'));
    if length(delims) < max(lfpcol, aligncol) - 1
        warning('lfp_makeCohFigs:shortline', ...
            'Skipping line %d, too short.', linenum);
        continue
    end
    delims = [0 delims length(line)+1];
    data(end+1).sessiondir = line(delims(sessioncol)+1 : delims(sessioncol+1)-1);
    data(end).cluster = line(delims(clustercol)+1 : delims(clustercol+1)-1);
    data(end).align = line(delims(aligncol)+1 : delims(aligncol+1)-1);
    data(end).lfp = line(delims(lfpcol)+1 : delims(lfpcol+1)-1);
end
fclose(fid);

%
% Make the coherograms
%
for fignum = 1:length(data)
    try
        lfp_declareGlobals;
        lfp_loadSetup;
        d = data(fignum);
        lfp_read2('preset', d.sessiondir,  ...
            {'Events.Nev', [LFPpre d.lfp '.' LFPext]} );
        lfp_AlignmentRef = str2num(d.align);
        % Assuming cluster names are in Rodent form (xxxn-TnCn, where xxx
        % is a string of 3 letters and n is a string of digits) or Naotaka
        % form (SEnCn).
        if monkeyflag
            Cidx = find(d.cluster == 'C');
            clustdir = d.sessiondir;
            lfp_add('preset', clustdir, ...
                {[d.cluster(1:Cidx-1) '.DWD']}, 'Naotaka Clusters (*.DWD)', false);
        else
            clustdir = fileparts(d.sessiondir);
            Tidx = find(d.cluster == 'T');
            Cidx = find(d.cluster(Tidx:end) == 'C') + Tidx - 1;
            trodenum = str2num(d.cluster(Tidx+1:Cidx-1));
            ext = sprintf('T%02i', trodenum);
            lfp_add('preset', clustdir, ...
                {[d.cluster(1:Tidx-2) '.' ext]}, ...
                'Rodent Clusters (*.Tnn, *.TTn)', false);
        end
        if isempty(evalstr)
            lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
            lfp_selectByDuration('and');
        else
            eval(evalstr);
        end
        lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
        clustnum = str2num(d.cluster(Cidx+1:end));
        lfp_createWave(@lfp_spike2wave, clustnum, ...
            'name', lfp_SpikeNames{clustnum});
        evtstr = [ num2str(lfp_AlignmentRef(1)) ];
        for k = 2:length(lfp_AlignmentRef)
            evtstr = [evtstr '-' num2str(lfp_AlignmentRef(k))];
        end
        pairstr = [lfp_SessionNames{1} '_' d.cluster '_' ...
            LFPpre d.lfp '_' evtstr '_'];
        for row = 1:size(nwk,1)
            hF = lfp_spec('coh', [], lfp_ActiveFilenums([1 2]), ...
                moving_win, 'avg', 'pad', padfactor, ...
                'logsig', 'rmdc', 'nw', nwk(row,1), 'k', nwk(row,2), ...
                'minspikes', 2*nwk(row,2));
            saveas(hF, sprintf('%s%g-%g.fig', ...
                pairstr, nwk(row,1), nwk(row,2) ));
            close(hF);
        end
    catch
        lfp_declareGlobals;
        if isempty(lfp_LogFileName)
            logname = 'lfp_lib.log';
            lfp_LogFileName = which(logname);
            if isempty(lfp_LogFileName)
                lfp_LogFileName = fullfile(pwd, logname);
            end
        end
        [msgstr, msgid] = lasterr;
        if exist('d') 
            if ~isempty(d)
                logmsg = sprintf('Error while processing %s clust %s align %s LFP %s\n%s\n%s', ...
                    d.sessiondir, d.cluster, d.align, d.lfp, ...
                    msgid, msgstr );
            else
                logmsg = sprintf('Error, empty d, %s\n%s', msgid, msgstr );
            end
        else
            logmsg = sprintf('Error %s\n%s', msgid, msgstr );
        end
        lfp_log(logmsg);
        disp(logmsg);
    end
end
