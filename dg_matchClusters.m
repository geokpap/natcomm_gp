function [clusterpairs, correl] = dg_matchClusters( ...
    fn1, fn2, thresh, varargin)
% Reads two files of cluster data from .MAT files <fn1>, <fn2>, and
% returns a list of clusters that match across the two files.
% Alternatively, numeric arrays containing the cluster data can be supplied
% instead of filenames for <fn1> and/or <fn2>, in which case the data
% supplied are used instead of reading data from a file.  The computation
% does not depend on how many wires worth of data the waveforms represent,
% so this function can be used with tetrode, stereotrode, or single
% electrode data.
%
% OUTPUT
% <clusterpairs> has two columns: the first contains a cluster ID from
% <fn1> and the second contains the last cluster ID from <fn2> whose
% waveform correlation coefficient with the cluster from <fn1> exceeds
% <thresh>.  Note that <thresh> should be adjusted so that there cannot be
% more than one pair with correl coeff > thresh.  <correl> is the matrix of
% correlation coefficients indexed by clusters from <fn1> on rows and
% clusters from <fn2> in columns.
%
% INPUT
% As for lfp_add 'Multiple Cut Cluster (*.MAT)', each file must contain an
% array with cluster ID in col 2 as first variable saved (does not care
% what its name is; col 1 is not used here but is tolerated for
% compatibility).  The remainder of the columns in the array must contain
% spike waveform sample data concatenated from all wires, i.e. if it is
% tetrode data with 32 samples per wire then there will 128 columns of
% waveform sample data.  It does not matter how many columns are actually
% supplied, as long as it is the same number of columns in each file.
%
% OPTIONS
% 'badwire', badwirelist - numbering the wires 1 - 4 in the same order in
%   which the sample channels appear in the columns of fn1 and fn2, all
%   sample values for any wires that appear in <badwirelist> are set to
%   zero in both cluster sets.  In the case of stereotrode data, address
%   the first wire as 1:2 and the second wire as 3:4.  WARNING: all other
%   things being equal, this will generally reduce the separation in
%   correlation values between repeated units and distinct units, to a
%   greater degree for each additional wire ignored.  However, it may
%   actually improve the separation in the situation for which it is
%   intended, i.e. where one wire has clearly gone bad.  Suffice it to say
%   that statistics and thresholds computed for four wires do not transfer
%   to three wires.
% 'ktfilenames' - passed through to dg_lookupSpikeWaveforms if necessary.
% 'plexonnames' - passed through to dg_lookupSpikeWaveforms if necessary.
% 'lookup' - this must be specified if fn1 and fn2 are ordinary 'Multiple
%   Cut Cluster (*.MAT)' files that do not contain waveform data.  The
%   waveform data is looked up from a matching Neuralynx tetrode data file
%   (if it can be found).
%
% EXAMPLES
% dg_matchClusters('file1', 'file2', .98, 'badwire', 4)
%   Matches clusters between file1 and file2 with a threshold of .98, but
%   ignoring any data on wire 4 in both files.  In this case, file1 and
%   file2 must contain full waveform data.
% dg_matchClusters('file1', 'file2', .99, 'lookup', 'badwire', [2 4])
%   Looks for the Neuralynx files containing the raw data that correspond
%   to file1 and file2, then compares clusters between them with a
%   threshold of .99 ignoring data on wire 2 AND wire 4 in both files.
%
% REFERENCE
% See p. 99 of Emondi AA, Rebrik SP, Kurgansky AV, Miller KD: Tracking
% neurons recorded from tetrodes across time.  J Neurosci Methods. 2004 May
% 30;135(1-2):95-105. 

%$Rev: 156 $
%$Date: 2012-09-18 18:42:20 -0400 (Tue, 18 Sep 2012) $
%$Author: dgibson $

lookupflag = false;
argnum = 1;
badwirelist = [];
dg_lsw_opts = {};
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'badwire'
            argnum = argnum + 1;
            badwirelist = varargin{argnum};
            if ~isa(badwirelist, 'numeric')
                error('dg_matchClusters:badoptval', ...
                    'The ''badwire'' option must be followed by a numeric array.');
            end
            dg_lsw_opts(end+1:end+2) = {'badwire', badwirelist};
        case 'ktfilenames'
            dg_lsw_opts(end+1) = {'ktfilenames'}; %#ok<*AGROW>
        case 'plexonnames'
            dg_lsw_opts(end+1) = {'plexonnames'};
        case 'lookup'
            lookupflag = true;
        otherwise
            error('dg_matchClusters:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Load the data into arrays s1 and s2:
if isa(fn1, 'numeric')
    s1 = fn1;
    fn1 = 'fn1 (numeric)';
else
    if lookupflag
        s1 = dg_lookupSpikeWaveforms(fn1, dg_lsw_opts{:});
    else
        S = load(fn1);
        fields = fieldnames(S);
        s1 = S.(fields{1});
        clear S;
    end
end
if isa(fn2, 'numeric')
    s2 = fn2;
    fn2 = 'fn2 (numeric)';
else
    if lookupflag
        s2 = dg_lookupSpikeWaveforms(fn2, dg_lsw_opts{:});
    else
        S = load(fn2);
        fields = fieldnames(S);
        s2 = S.(fields{1});
        clear S;
    end
end

if size(s1,2) ~= size(s2,2)
    error('dg_matchClusters:baddata', ...
        '<fn1> and <fn2> must contain the same number of data columns.' );
end

if ~isempty(badwirelist)
    if mod(size(s1,2) - 2, 4) ~= 0
        error('dg_matchClusters:baddata2', ...
        'The total number of sample points must be divisible by 4.' );
    end
    numptsperwire = round((size(s1,2) - 2) / 4);
    startingbadcols = numptsperwire * (badwirelist - 1) + 3;
    endingbadcols = numptsperwire * badwirelist + 2;
    badcols = [];
    for k = 1:numel(badwirelist)
        badcols = [badcols startingbadcols(k):endingbadcols(k)];
    end
    s1(:, badcols) = 0;
    s2(:, badcols) = 0;
end

waves1 = avgClustWaves(s1);
waves2 = avgClustWaves(s2);
figure;
subplot(2, 1, 1);
plot(waves1(:, 2:end)');
for k = 1:size(waves1, 1)
    legendstr{k} = num2str(waves1(k,1));
end
legend(legendstr);
title(sprintf('%s', fn1), 'Interpreter', 'none');
legendstr={};
subplot(2, 1, 2);
plot(waves2(:, 2:end)');
for k = 1:size(waves2, 1)
    legendstr{k} = num2str(waves2(k,1));
end
legend(legendstr);
title(sprintf('%s', fn2), 'Interpreter', 'none');
magnitudes1 = sqrt(sum(waves1(:, 2:end) .^ 2, 2));
magnitudes2 = sqrt(sum(waves2(:, 2:end) .^ 2, 2));
clusterpairs = zeros(size(waves1,1), 2);
for waverow1 = 1:size(waves1,1)
    for waverow2 = 1:size(waves2,1)
        correl(waverow1, waverow2) = ...
            sum(waves1(waverow1, 2:end) .* waves2(waverow2, 2:end)) ...
            / (magnitudes1(waverow1) * magnitudes2(waverow2));
        if abs(correl(waverow1, waverow2)) > thresh
            clusterpairs(waverow1, :) = ...
                [waves1(waverow1,1) waves2(waverow2,1)];
        end
    end
end
clusterpairs(clusterpairs(:,1) == 0, :) = [];


function avgwaves = avgClustWaves(s)
% Returns cluster number in col 1 and averaged sample data in remaining
% columns.
avgwaves = [];
clusters = unique(s(:,2))';
for clustnum = clusters
    clustselect = s(:,2) == clustnum;
    avgwaves(end+1, :) = [ clustnum mean(s(clustselect,3:end)) ];
end

