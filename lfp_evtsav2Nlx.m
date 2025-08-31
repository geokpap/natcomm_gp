function lfp_evtsav2Nlx(sessiondir)
%lfp_evtsav2Nlx(sessiondir)
% For each *.evtsav file in sessiondir, creates a *.nev file containing the
% same events.  Does not include "trial parameters", but only those event
% markers that eventually end up in lfp_Events.  The output filename is
% composed from the original *.evtsav filename such that it does not
% overwrite any existing file.

%$Rev: 333 $
%$Date: 2014-10-16 18:56:27 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

evtstrfmt = 'RecID: 4098 Port: 0 TTL Value: 0x%04X';
evtstrlen = length(sprintf(evtstrfmt, 0));
hdr = {
    '######## Ceci n''est pas un Neuralynx Data File Header'
    '## created by lfp_evtsav2Nlx'
    };
filenames = dir(sessiondir);
myfilenames = {};
for k = 1:length(filenames)
    if ismember(filenames(k).name, {'.', '..'}) ...
            || (filenames(k).bytes == 0)
        continue
    end
    [p, name, ext] = fileparts(filenames(k).name);
    if strcmpi(ext, '.evtsav')
        load('-mat', fullfile(sessiondir, filenames(k).name));
        TTLids = lfp_save_events(:,2) + 32768;
        TSs = round(lfp_save_events(:,1) * 1e+6);
        evtstrings = char(zeros(size(TTLids,1), evtstrlen));
        thing2 = 119 * ones(1, size(TTLids,1));
        thing4 = zeros(8, size(TTLids,1));
        for k = 1:size(TTLids,1)
            evtstrings(k,:) = sprintf(evtstrfmt, TTLids(k));
        end
        vnum = 0;
        outfname = fullfile(sessiondir, [name '.nev']);
        while exist(outfname)
            vnum = vnum+1;
            outfname = fullfile(sessiondir, ...
                sprintf('%s_%d.nev', name, vnum));
        end
        hdr{3} = sprintf('## Time Opened: %s', datestr(now));
        if ~ispc
            error('lfp_evtsav2Nlx:notPC', ...
                'Neuralynx format files can only be created on Windows machines.');
        end
        Mat2NlxEV(outfname, ...
            TSs', thing2, TTLids', thing4, cellstr(evtstrings), size(TTLids,1));
    end
end
