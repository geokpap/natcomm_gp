function dg_localAvgRefIfNeeded(srcroot, destroot, sessionID, ...
    trodegroups, varargin)
%dg_LocalAvgRefIfNeeded(srcroots, destroot, sessionID)
% Check to see if there is a complete set of local average referenced files
% in <sessionID> under directory <destroot>.  If there is not a complete
% set, generate an entire new set, overwriting any pre-existing files with
% the same names.
%   The downsampling factor is extracted from the filenames of the form
% csc<CH>_down<N>.MAT, where <CH> is the channel number and <N> is the
% downsampling factor.  If there is a 'dg_downsampleIfNeeded.mat' file in
% <srcroot>, then the downsampling factor is read from that.  Otherwise,
% the same rules used by dg_downsampleIfNeeded are used to determine it.
%INPUTS
% srcroot: directory containing session directories containing CSC files
%   referenced as recorded.
% destroot: directory containing session directories containing local
%   average referenced CSC files and the reference files themselves.
% sessionID: exactly as listed in <srcroot>.
% trodegroups: a cell array containing one numerical array of electrode
%   numbers in each cell.  Each cell defines one local average reference
%   group.
%OUTPUTS
% All output is directed to the following set of files:
%   <prefix><N>.mat files created by dg_makeLocalAvgRefs in directory
%       fullfile(destroot, sessionID).
%   <infilename>-<refname>.mat files created by dg_subtractLocalAvgRefs
%       in directory fullfile(destroot, sessionID).
%   File named <prefix>_groups.mat containing <trodegroups>.  There is
%       exactly one full set of local avg files for one value of <prefix>;
%       if less than a full set is found, any old files are overwritten.
%       Default value of <prefix> is 'ref'.
%OPTIONS
% 'lfpprefix', lfpprefix - when constructing LFP filenames, replaces the
%   'csc' in default formats 'csc%d.mat' and 'csc%d.ncs' with <lfpprefix>.
% 'mintrodes', N - sets the minimum number electrodes required to process a
%   group.  Groups that contain fewer than <N> channels are simply skipped.
%   Default N = 3;
% 'prefix', prefix - Sets value of <prefix> submitted to
%   dg_makeLocalAvgRefs and dg_subtractLocalAvgRefs.

%$Rev: 215 $
%$Date: 2015-03-27 01:05:37 -0400 (Fri, 27 Mar 2015) $
%$Author: dgibson $

lfpprefix = 'csc';
prefix = 'ref';
mintrodes = 3; 
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'lfpprefix'
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('dg_localAvgRefIfNeeded:lfpprefix', ...
                    '''lfpprefix'' option requires a value.');
            end
            lfpprefix = varargin{argnum};
            if ~ischar(lfpprefix)
                error('dg_localAvgRefIfNeeded:lfpprefix2', ...
                    '''lfpprefix'' option requires a string value.');
            end
        case 'mintrodes'
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('dg_localAvgRefIfNeeded:mintrodes', ...
                    '''mintrodes'' option requires a value.');
            end
            mintrodes = varargin{argnum};
            if ~isnumeric(mintrodes)
                error('dg_localAvgRefIfNeeded:mintrodes2', ...
                    '''mintrodes'' option requires a numeric value.');
            end
        case 'prefix'
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('dg_localAvgRefIfNeeded:prefix', ...
                    '''prefix'' option requires a value.');
            end
            prefix = varargin{argnum};
            if ~ischar(prefix)
                error('dg_localAvgRefIfNeeded:prefix2', ...
                    '''prefix'' option requires a string value.');
            end
        otherwise
            error('dg_localAvgRefIfNeeded:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

requiredtrodes = [];
for grpidx = 1:length(trodegroups)
    if numel(trodegroups{grpidx}) >= mintrodes
        requiredtrodes = [
            requiredtrodes
            reshape(trodegroups{grpidx}, [], 1)
            ]; %#ok<AGROW>
    else
        trodegroups{grpidx} = [];
    end
end
requiredtrodes = unique(requiredtrodes);
sessiondir = fullfile(srcroot, sessionID);
N = dg_findDownsamplingFactor(sessiondir, requiredtrodes);

% check to see if there is a complete set of local avg referenced files
outputfiles = {};
for grpidx = 1:length(trodegroups)
    for trodeidx = 1:length(trodegroups{grpidx})
        if isnan(N)
            outputfiles{end+1,1} = sprintf( 'csc%d-%s%d', ...
                trodegroups{grpidx}(trodeidx), prefix, grpidx ); %#ok<AGROW>
        else
            % The '.*' is a wildcard to match the optional downsampling
            % factor that may precede <prefix>:
            outputfiles{end+1,1} = sprintf( ...
                'csc%d_down%d-.*%s%d', ...
                trodegroups{grpidx}(trodeidx), N, prefix, grpidx ); %#ok<AGROW>
        end
    end
end
destdir = fullfile(destroot, sessionID);
gottem = true;
for fileidx = 1:length(outputfiles)
    if isempty(dg_fileExists(destdir, outputfiles{fileidx}))
        gottem = false;
        break
    end
end

if gottem
    fprintf('dg_LocalAvgRefIfNeeded: all files exist for %s\n', destdir);
    return
else
    % construct <localavgrefs> from <trodegroups>
    fprintf('dg_LocalAvgRefIfNeeded: creating files for %s\n', destdir);
    localavgrefs = cell(size(trodegroups));
    for grpidx = 1:length(trodegroups)
        for trodeidx = 1:length(trodegroups{grpidx})
            if isnan(N)
                % These are raw-sampling-rate files, so they could be either
                % .mat or .ncs
                filename = sprintf('%s%d.mat', lfpprefix, ...
                    trodegroups{grpidx}(trodeidx));
                result = dg_fileExists(sessiondir, filename);
                if isempty(result)
                    filename = sprintf('%s%d.ncs', lfpprefix, ...
                        trodegroups{grpidx}(trodeidx));
                    if ~exist(fullfile(sessiondir, filename), 'file')
                        error('dg_LocalAvgRefIfNeeded:missing', ...
                            'Required file %s is missing', ...
                            fullfile(sessiondir, filename));
                    end
                else
                    if length(result) > 1
                        warning('dg_LocalAvgRefIfNeeded:multi', ...
                            'Multiple versions of filename: %s', ...
                            dg_thing2str(result));
                    end
                    filename = result{1};
                end
            else
                filename = sprintf('csc%d_down%d\\.mat', ...
                    trodegroups{grpidx}(trodeidx), N);
                result = dg_fileExists(sessiondir, filename);
                if isempty(result)
                    error('dg_LocalAvgRefIfNeeded:missing2', ...
                        'Required file %s is missing', ...
                        fullfile(sessiondir, filename));
                else
                    if length(result) > 1
                        warning('dg_LocalAvgRefIfNeeded:multi2', ...
                            'Multiple versions of filename: %s', ...
                            dg_thing2str(result));
                    end
                    filename = result{1};
                end
            end
            localavgrefs{grpidx}{trodeidx} = filename;
        end
    end
    
    save(fullfile(destdir, sprintf('%s_groups.mat', prefix)), ...
        'trodegroups');
    
    % Make reference files
    dg_makeLocalAvgRefs(localavgrefs, sessiondir, prefix);
    
    % Make referenced files
    dg_subtractLocalAvgRefs(localavgrefs, sessiondir, prefix);
end

