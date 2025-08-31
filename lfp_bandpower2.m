function [result, t] = lfp_bandpower2(samples, moving_win, freqlims, intrp)
%result = lfp_bandpower2(samples, moving_win, freqlims, intrp)
%	Modified version of lfp_bandpower for more general utility use.
%   Accepts literal sample data instead of a filenum; <samples> must be
%   a column vector.
%   The call to lfp_interpolateWave constrains this func only to work on
%   <samples> of same length as an entire CSC channel; consequently this
%   amount of storage is pre-allocated before calling lfp_interpolateWave.
%   Accepts multiple bands on multiple rows of <freqlims>, returns
%   corresponding multiple results in multiple columns of <result>.
%   <intrp> is a logical value that controls whether or not the result gets
%   interpolated by lfp_interpolateWave.  <t> is a vector of timestamps of
%   the centers of the moving windows.
%   
%   NOTE: Cannot disable action of lfp_FreqLim due to call to
%   lfp_stspecgram, which respects lfp_FreqLim.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if size(samples,2) ~= 1
    error('lfp_bandpower2:badsamples', ...
        '<samples> must be a column vector.' );
end
if size(samples,1) ~= length(lfp_TimeStamps) * lfp_SamplesPerFrame
    error('lfp_bandpower2:badsamples2', ...
        '<samples> must be same length as CSC channels.' );
end

approx_f = 0 : 1/moving_win(1) : 1/(2*lfp_SamplePeriod);
for bandnum = 1:size(freqlims,1)
    selfreqs = find(approx_f >= freqlims(bandnum, 1) ...
        & approx_f <= freqlims(bandnum, 2));
    if ~isempty(lfp_FreqLim)
        if any(approx_f(selfreqs) < lfp_FreqLim(1)) || ...
                any(approx_f(selfreqs) > lfp_FreqLim(2))
            warning('lfp_bandpower2:trimmingfreqs', ...
                'Eliminating specified freqs for band %d that are outside lfp_FreqLim', ...
                bandnum );
        end
        selfreqs = find(approx_f(selfreqs) >= lfp_FreqLim(1) ...
            & approx_f(selfreqs) <= lfp_FreqLim(2));
    end
    if isempty(selfreqs)
        error('lfp_bandpower2:nofreqs', ...
            'No frequencies selected for band %d, check <freqlims> and lfp_FreqLim.', ...
            bandnum );
    end
end

approx_num_windows = length(samples) ...
    * lfp_SamplePeriod / moving_win(2);
if approx_num_windows > 10000
    warning('lfp_bandpower2:manywindows', ...
        'You are computing approximately %d points, at about 1 minute per 10000', ...
        approx_num_windows);
end

[S f wsl] = lfp_stspecgram(@hamming, ...
    samples, moving_win, ...
    1/lfp_SamplePeriod, [], 0, 0);

for bandnum = 1:size(freqlims,1)
    selfreqs = find(f >= freqlims(bandnum, 1) & f <= freqlims(bandnum, 2));
    disp('The frequencies to be summed are:');
    disp(reshape(f(selfreqs),[],1));
    wavedata(:,bandnum) = sum(S(:, selfreqs), 2);
end
t = zeros(size(wsl));
% Find timestamps of start of each window:
hWaitBar = waitbar(0, '', 'Name', 'Finding timestamps of windows');
for idx = 1:length(wsl)
    waitbar(idx/length(wsl), hWaitBar );
    t(idx) = lfp_index2time(wsl(idx));
end
close(hWaitBar);
% Convert to middle of each window:
t = t + moving_win(1)/2;

result = zeros(numel(samples), size(freqlims,1));
for bandnum = 1:size(freqlims,1)
    if intrp
        result(:,bandnum) = lfp_interpolateWave(t, wavedata(:,bandnum));
    else
        result = wavedata;
    end
end

