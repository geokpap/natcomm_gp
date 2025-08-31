function dg_spikes2nev(outputfilename, spikefiles, nevfiles)
%dg_spikes2nev(outputfilename, spikefiles, nevfiles)
% Merges events recorded in spike files with events recorded in events
% files and saves the result as a Neuralynx events file.
%INPUT
% outputfilename: pathname to output file; may be absolute or relative.
% spikefiles: a two-column cell array; first column contains TTL codes for
%   the events recorded in the spike files, the second column contains
%   pathnames to the spike files, so that <spikefiles{k,1}> is the TTL code
%   to use for the events in <spikefiles{k,2}>.  <spikefiles{k,2}> may be
%   absolute or relative, but if it's relative, it is NOT read from the
%   current working directory, but rather from the same directory as
%   <nevfiles{1}>.  Note that this function does not do anything to the
%   strobe bit, so if you want TTLs that have the strobe bit set, it is
%   necessary to add 2^15 from your desired value.  For example:
%       spikefiles = {
%           17 + 2^15   'spikes1.nse'
%           31 + 2^15   'spikes2.nse' };
%   The spike files are assumed to be single-channel spike files, as in the
%   example.  WARNING: do not ask how it is possible that the conversion
%   from unstrobed to strobed is accomplished by the same mathematical
%   operation as the conversion from strobed to unstrobed.  It is a mystery
%   that passeth understanding and should be avoided like uttering the
%   tetragrammaton.
% nevfiles: a cell string array of Neuralynx event files to merge with the
%   events in the spike files.

%$Rev: 134 $
%$Date: 2011-11-17 16:25:41 -0500 (Thu, 17 Nov 2011) $
%$Author: dgibson $

if ~iscell(spikefiles) || ~iscell(nevfiles)
    error('dg_spikes2nev:argtype1', ...
        'Both <spikefiles> and <nevfiles> must be cell arrays.');
end
if ~ischar(outputfilename)
    error('dg_spikes2nev:argtype2', ...
        '<outputfilename> must be a string.');
end

events = zeros(0, 2);
evtstr = cell(0, 1);
header = '';
for fileidx = 1:length(nevfiles)
    try
        [TS, TTL, ES, Hdr] = dg_readEvents(nevfiles{fileidx});
    catch e
        logmsg = sprintf('%s\n%s', ...
            e.identifier, e.message);
        for stackframe = 1:length(e.stack)
            logmsg = sprintf('%s\n%s\nline %d', ...
                logmsg, e.stack(stackframe).file, e.stack(stackframe).line);
        end
        warning('dg_spikes2nev:badnev', ...
            'Could not read events file %s\n%s', nevfiles{fileidx}, logmsg);
        TS = [];
        TTL = [];
        ES = {};
        Hdr = {'######## Neuralynx Data File Header'};
    end
    if fileidx == 1
        [pathstr, name, ext] = fileparts(nevfiles{fileidx}); %#ok<*NASGU,*ASGLU>
        if isempty(pathstr)
            srcdir = pwd;
        else
            srcdir = pathstr;
        end
        header = Hdr;
    end
    events = [ events
        reshape(TS, [], 1) reshape(TTL, [], 1) ]; %#ok<*AGROW>
    evtstr = [ evtstr
        reshape(ES, [], 1) ];
end
for fileidx = 1:size(spikefiles,1)
    [pathstr, name, ext] = fileparts(spikefiles{fileidx,2});
    if isempty(pathstr)
        filename = fullfile(srcdir, spikefiles{fileidx,2});
    else
        filename = spikefiles{fileidx,2};
    end
    TS = Nlx2MatSpike_411(filename, [1 0 0 0 0], 0, 1);
    events = [ events
        reshape(TS, [], 1) repmat(spikefiles{fileidx,1}, length(TS), 1) ];
    ES = repmat(sprintf('dg_spikes2nev event 0x%04X', ...
            spikefiles{fileidx,1}), length(TS), 1);
    evtstr = [ evtstr
        cellstr(ES) ];
end
[events, idx] = sortrows(events);
evtstr = evtstr(idx);
Hdr = header(1);
Hdr{2,1} = sprintf( ...
    '## File Name: %s ', outputfilename);
Hdr{3,1} = sprintf( ...
    '## srcdir: %s ', srcdir);
infilenames = [ reshape(nevfiles, [], 1)
    spikefiles(:,2) ];
for filenum = 1:length(infilenames)
    Hdr{end+1} = sprintf('## %s ', infilenames{filenum});
end
dg_writeNlxEvents(outputfilename, events(:,1), events(:,2), evtstr, Hdr);


