function [sampidx, samples, params] = dg_kilosort2snapshots(sessiondir, ...
    numpresamples, numpostsamples)
% Extracts multichannel spike snapshots in analogous fashion to a Neuralynx
% spike file.
%INPUTS
% sessiondir: path to directory containing a collection of 'kilosort'
%   output files.
% numpresamples: the number of samples to return before the nominal spike
%   time.
% numpostsamples: the number of samples to return after the nominal spike
%   time.
%OUTPUTS
% sampidx:  a cell row vector containing one cell for each spike cluster.
%   Each cell contains a column vector of the sample number at which each
%   spike in the cluster nominally occured.
% samples: a cell row vector containing one cell for each spike cluster.
%   Each cell contains the multichannel spike waveforms for all spikes in
%   the cluster in (channels, samples, spikes) format.
% params: a cell string vector containing the contents of 'params.py'.
%NOTES
%   The total number of samples returned per channel for each spike is
% <numsamp> = <numpresamples> + <numpostsamples> + 1.
%   Kilosort can separate (partially) superimposed spikes; it is therefore
% important when analyzing individual spike shapes to look for and ignore
% any spikes that occur inside another spike's time window.

%$Rev: 259 $
%$Date: 2017-08-02 14:37:57 -0400 (Wed, 02 Aug 2017) $
%$Author: dgibson $

sampidx = []; %#ok<NASGU>
samples = []; %#ok<NASGU>

params = fileread(fullfile(sessiondir, 'params.py'));
% We trust that 'params.py' is exactly as specified, and simply evaluate
% its contents line by line to get the values of the variables dat_path,
% n_channels_dat, dtype, offset, sample_rate, and hp_filtered.
params = strsplit(params, '\n');
for k = 1:length(params)
    eval(lower(params{k}));
end
sampidx = cell(1, n_channels_dat);
samples = cell(1, n_channels_dat);

spike_times = double(readNPY(fullfile(sessiondir, 'spike_times.npy')));
spike_templates = readNPY(fullfile(sessiondir, 'spike_templates.npy'));
clustIDs = unique(spike_templates);
numclusts = length(clustIDs);
if ~isequal(clustIDs, 1 : numclusts)
    warning('dg_kilosort2snapshots:clustIDs\n%s', ...
        'cluster IDs are not consecutive integers:', ...
        dg_canonicalSeries(clustIDs));
end

fprintf('Reading CSC data...\n');
rawdatafile = fullfile(sessiondir, dat_path);
fid = fopen(rawdatafile);
if fid == -1
    error('dg_kilosort2snapshots:rawdatafile', ...
        'Could not open %s.', rawdatafile);
else
    try
        CSCdata = fread(fid, dtype);
    catch e
        fclose(fid);
        fprintf( ...
            'error while reading %s', rawdatafile);
        rethrow(e);
    end
    if ismember(fid, fopen('all'))
        fclose(fid);
    end
end
fprintf('Finished reading CSC data.\n');
% CSCdata(channelidx, sampleidx)
CSCdata = reshape(CSCdata, n_channels_dat, []);

% Sort spike times and waveforms by cluster
relsampidx = -numpresamples : numpostsamples;
for clustidx = 1:numclusts
    sampidx{clustidx} = reshape( ...
        spike_times(spike_templates == clustIDs(clustidx)), ...
        [], 1 );
    % Retrieve waveforms in (channels, samples, spikes) format. (This just
    % gets too convoluted when I think about parallelizing it.)
    samples{clustidx} = NaN( n_channels_dat, length(relsampidx), ...
        length(sampidx{clustidx}) );
    for spkidx = 1:length(sampidx{clustidx})
        absampidx = relsampidx + sampidx{clustidx}(spkidx);
        samples{clustidx}(:, :, spkidx) = CSCdata(:, absampidx);
    end
end


