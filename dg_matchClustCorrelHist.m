function correlations = dg_matchClustCorrelHist(varargin)
%correlations = dg_matchClustFalseHist(sessiondir1)
%correlations = dg_matchClustFalseHist(sessiondir1, sessiondir2)
%correlations = dg_matchClustFalseHist('ktfilenames')
%
% If sessiondir2 is not given, finds all Multiple Cut Cluster .MAT files
% from the session represented by sessiondir1 (i.e. files in the parent
% directory of sessiondir1 with filenames of the form 'sessionid-Tnn.MAT',
% e.g. 'ACQ06-T03.MAT'), and computes the waveform correlations between all
% pairs of clusters from different files, i.e. between all independent
% pairs of distinct units. (Clusters from the same electrode could
% potentially represent the same unit, and necessarily have anticorrelated
% waveforms because of the disjoint nature of cluster cutting.) 
%
% If sessiondir2 is given, performs the same between all pairs of clusters
% from different sessions in files of the same name, i.e. between pairs of
% clusters recorded on the same electrode in different sessions.
%
% In either case, displays a histogram of the distribution of correlation
% values, and returns the complete set of values in <correlations>.
%
%OPTIONS
% In this perverted function, the options can go anywhere in the arg list.
% 'ktfilenames' - filenames are of the form 'sessionid_nn*.MAT', instead of
%   'sessionid-Tnn.MAT', where <nn> is a 2-digit electrode number that may
%   be followed by any number of characters of gobbledygook.  A desperate
%   search to find the preferred set of cluster files in the subdirectories
%   of the session dir is also undertaken, rather than assuming they are in
%   the animal directory.

%$Rev: 32 $
%$Date: 2009-05-26 19:11:53 -0400 (Tue, 26 May 2009) $
%$Author: dgibson $

argnum = 1;
args2delete = [];
ktfflag = false;
myoptions = {};
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'badwire'
            args2delete(end+1) = argnum;
            argnum = argnum + 1;
            badwirelist = varargin{argnum};
            args2delete(end+1) = argnum;
            myoptions(end+1:end+2) = {'badwire', badwirelist};
        case 'ktfilenames'
            ktfflag = true;
            args2delete(end+1) = argnum;
            myoptions(end+1) = varargin(argnum);
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];

numsessions = length(varargin);
if ~ismember(numsessions, [1 2])
    error('dg_matchClustCorrelHist:numsessions', ...
        'You may not specify more than two sessions for comparison.');
end

% <clustfileregexp> matches the part after the sessionID
if ktfflag
    clustfileregexp = '_\d\d.*\.MAT';
else
    clustfileregexp = '-T\d\d+\.MAT';
end

% Find spike files, get waveforms
for snum = 1:numsessions
    sessiondir = varargin{snum};
    if ktfflag
        if exist([sessiondir '\onlytrialcut\ccw'], 'dir')
            clustfiledir = [sessiondir '\onlytrialcut\ccw'];
        elseif exist([sessiondir '\onlytrialcut'], 'dir')
            clustfiledir = [sessiondir '\onlytrialcut'];
        else
            error('Burma!');
        end
        MATfiles = dir(fullfile(clustfiledir, '*.MAT'));
        inds = regexp(sessiondir, '\');
        animaldir = sessiondir(1:inds(end-1)-1);
        sessionID = sessiondir(regexpi(sessiondir, ...
            'acq\d\d'):regexpi(sessiondir, 'acq\d\d')+4);
    else
        [animaldir, sessionID] = fileparts(sessiondir);
        MATfiles = dir(fullfile(animaldir, [sessionID '-T*.MAT']));
        clustfiledir = animaldir;
    end
    % Eliminate any files that don't match the convention for Multiple Cut
    % Cluster .MAT filenames:
    mismatch = [];
    for k = 1:length(MATfiles)
        if isempty(regexpi(MATfiles(k).name, [sessionID clustfileregexp]))
            mismatch(end+1) = k;
        end
    end
    MATfiles(mismatch) = [];
    % Retrieve the spike waveforms
    spikewaves{snum} = cell(size(MATfiles));
    for k = 1:length(MATfiles)
        delim1 = regexpi(MATfiles(k).name, clustfileregexp);
        if ktfflag
            trodestr = MATfiles(k).name(delim1+1:delim1+2);
        else
            sessionID = MATfiles(k).name(1 : delim1 - 1);
            trodestr = MATfiles(k).name(delim1 + 2 : end);
            % chop any non-digits from end:
            delim2 = regexp(trodestr, '\D');
            if ~isempty(delim2)
                trodestr = trodestr(1:delim2(1));
            end
        end
        trodenums{snum}(k) = str2num(trodestr);
        spikewaves{snum}{k} = dg_lookupSpikeWaveforms( ...
            fullfile(clustfiledir, MATfiles(k).name), myoptions{:} );
    end
end

correlations = [];
thresh = 0.9;
if numsessions == 1
    for k = 1 : (length(spikewaves{1}) - 1)
        for m = (k+1) : length(spikewaves{1})
            [clusterpairs, correl] = dg_matchClusters( ...
                spikewaves{1}{k}, spikewaves{1}{m}, thresh, myoptions{:} );
            correlations = [correlations; correl(:)];
        end
    end
else % numsessions = 2
    for k = 1 : length(spikewaves{1})
        if length(trodenums) < 2 || isempty(trodenums{2})
            m = [];
        else
            m = find(ismember(trodenums{2}, trodenums{1}(k)));
        end
        if isempty(m)
            warning('dg_matchClustCorrelHist:nomatch', ...
                'There is no electrode %d file in %s, skipping', ...
                trodenums{1}(k), varargin{2} );
        else
            [clusterpairs, correl] = dg_matchClusters( ...
                spikewaves{1}{k}, spikewaves{2}{m}, thresh, myoptions{:} );
            correlations = [correlations; correl(:)];
        end
    end
end

figure;
bins=(0:200)/200;
hist(correlations, bins);
