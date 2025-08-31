function dg_computeCSCStats(sessiondir, varargin)
%INPUT
% sessiondir: an absolute or relative pathname for a directory containing
%   raw Neuralynx session data files.
%OUTPUT
% Saves output to a file in <sessiondir> named 'computeCSCStats.mat',
% containing the following variables.  If there is already a
% 'computeCSCStats.mat', the old one is moved to 'computeCSCStats_bak.mat'
% before creating the new one.  However, this is done after analyzing each
% CSC file, so if everything runs normally to completion, then
% 'computeCSCStats.mat' contains the stats for all the files, and
% 'computeCSCStats_bak.mat' contains the exact same stats for all except
% the last file, so you should not rely on this feature to save a
% pre-existing 'computeCSCStats.mat'.
%   files: a cell array of the names of the files that were analyzed.
%   CSCStats: a vector with one element per file listed in <files>, each
%     element of which is a structure containing the following fields, all
%     of which are numeric:
%       ADmaxval: the -ADMaxValue stored in the corresponding file's
%         header.
%       bitvolts: the -ADBitVolts from the file header.
%       disorderedTSidx: indices (i.e. record numbers) of the first record
%         in any pair where the timestamp of the second record is earlier
%         than the timestamp of the first.
%       max: the absolute maximum sample value in the file.
%       min: the absolute minimum sample value in the file.
%       dc: the DC offset in volts, computed as the average value of the
%           entire file.
%       medianpeak: the median absolute value of the peaks and valleys in
%           the AC waveform.
%       fractionClipped: the fraction of samples that are equal to the
%           maximum or minimum possible AD values.
%     If any of the above fields contains NaN, that means the information
%     could not be extracted from the file.
%OPTIONS
% 'destdir', destdir - saves output file to <destdir> instead of to
%   <sessiondir>.

destdir = '';
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'destdir'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        otherwise
            error('dg_computeCSCStats:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(destdir)
    destdir = sessiondir;
end

files = dir(fullfile(sessiondir, '*.ncs'));
FILES = dir(fullfile(sessiondir, '*.NCS'));
files = [files
    FILES];

for fileidx = 1:length(files)
    fprintf('Reading %s...\n', files(fileidx).name);
    [TS, Samples, Hdr] = dg_readCSC( ...
        fullfile(sessiondir, files(fileidx).name) );
    ADMaxValMatch = regexpi(Hdr, '-ADMaxValue\s+(\d+)', 'tokens');
    istheone = ~cellfun(@isempty, ADMaxValMatch);
    if ~any(istheone)
        CSCStats(fileidx).ADmaxval = NaN;
    else
        ADMaxValStr = ADMaxValMatch{istheone}{1};
        CSCStats(fileidx).ADmaxval = str2double(ADMaxValStr);
    end
    ADBitVMatch = regexpi(Hdr, '-ADBitVolts\s+([.\d]+)', 'tokens');
    istheone = ~cellfun(@isempty, ADBitVMatch);
    if ~any(istheone)
        CSCStats(fileidx).bitvolts = NaN;
    else
        ADBitVStr = ADBitVMatch{istheone}{1};
        CSCStats(fileidx).bitvolts = str2double(ADBitVStr);
    end
    CSCStats(fileidx).disorderedTSidx = find(diff(TS) < 0); %#ok<*AGROW>
    CSCStats(fileidx).max = max(Samples(:));
    CSCStats(fileidx).min = min(Samples(:));
    if CSCStats(fileidx).max == CSCStats(fileidx).ADmaxval || ...
            CSCStats(fileidx).min == -CSCStats(fileidx).ADmaxval - 1
        numclipped = sum( Samples(:) == CSCStats(fileidx).ADmaxval | ...
            Samples(:) == -CSCStats(fileidx).ADmaxval - 1 );
        CSCStats(fileidx).fractionClipped = numclipped / numel(Samples);
    else
        CSCStats(fileidx).fractionClipped = 0;
    end
    if isnan(CSCStats(fileidx).bitvolts)
        warning('computeCSCStats:dc', ...
            'Cannot evaluate DC because there is no bitvolts value');
    end
    CSCStats(fileidx).dc = CSCStats(fileidx).bitvolts * mean(Samples(:));
    pospeakidx = dg_findpks(Samples(:));
    negpeakidx = dg_findpks(-Samples(:));
    CSCStats(fileidx).medianpeak = median(abs( ...
        Samples([pospeakidx negpeakidx]) - CSCStats(fileidx).dc ));
    if exist(fullfile(sessiondir, 'computeCSCStats.mat'), 'file')
        [status, message] = movefile( ...
            fullfile(sessiondir, 'computeCSCStats.mat'), ...
            fullfile(sessiondir, 'computeCSCStats_bak.mat') );
        if ~status
            warning('computeCSCStats:move', ...
                'Could not move output file to backup:\n%s', ...
                message );
        end
    end
    fprintf('Saving stats through %s...\n', files(fileidx).name);
    save(fullfile(destdir, 'computeCSCStats.mat'), 'files', 'CSCStats');
end
