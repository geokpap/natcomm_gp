function [result, units] = lfp_nhtfreq(filenum, varargin)
%lfp_nhtfreq(filenum) for use with lfp_createWave; normalized Hilbert
% transform calculation of instantaneous frequency from the waveform in
% <filenum>.  Finds runs of consecutive non-NaN and processes each such run
% separately.
%OPTIONS
% 'verbose'

%$Rev: 197 $
%$Date: 2011-01-20 00:37:59 -0500 (Thu, 20 Jan 2011) $
%$Author: dgibson $

global lfp_Samples lfp_SamplePeriod

verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_nhtfreq:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

result = NaN(size(lfp_Samples{filenum}));
units = 'Hz';

waveidx = reshape(~isnan(lfp_Samples{filenum}), 1, []);
wavestarts = find([ waveidx(1) waveidx(2:end) & ~waveidx(1:end-1) ]);
waveends = find([ waveidx(1:end-1) & ~waveidx(2:end) waveidx(end) ]);
if length(wavestarts) ~= length(waveends)
    error('lfp_nhtfreq:botch', ...
        'Dan botched something again :-(');
end
for wavenum = 1:length(wavestarts)
    if verboseflag
        fprintf('starting #%d of %d wave segments\n', ...
            wavenum, length(wavestarts) );
    end
    startsample = wavestarts(wavenum);
    endsample = waveends(wavenum);
    result(startsample:endsample) = dg_nht(lfp_Samples{filenum}( ...
        startsample:endsample ), 'pchip') / (lfp_SamplePeriod*2*pi);
end
