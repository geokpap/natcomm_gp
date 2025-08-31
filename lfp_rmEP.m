function [wave, units]  = lfp_rmEP(filenum, xstart, y, varargin)
%wave = lfp_rmEP(filenum, xstart, y)

% Function intended for use with lfp_createWave. Subtracts the Evoked
% Potential specified by <y> from the wave in <filenum> and returns the
% result.  <y> is reshaped as a vector regardless of its original shape.
% <xstart> represents the time of the first sample in <y> relative to
% lfp_AlignmentRef.  Any time points in <wave> that are not within the time
% interval covered by the EP for some trial are given the value NaN.  If a
% time point is within the EP interval for more than one trial, a warning
% is issued, and the EP-removed trace of the earlier of the two trials gets
% truncated at the point where the trace for the later trial starts.  If a 
% trial does not contain an lfp_AlignmentRef event, it is
% skipped with a warning.
%OPTIONS
% 'taper' - applies a Hanning window (excluding the zero-weighted end
%   samples) to <y> before subtracting it, and fills the intervening
%   intervals with an unmodified copy of the waveform in <filenum> instead
%   of NaN.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~ismember(filenum, lfp_ActiveFilenums)
    error('lfp_rmEP:badclust', ...
        'There is no file number %d', filenum);
end
if numel(xstart) ~= 1
    error('lfp_rmEP:baddata', ...
        '<xstart> must be a scalar');
end

taperflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'taper'
            taperflag = true;
        otherwise
            error('lfp_rmEP:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

y = reshape(y,1,[]);
if taperflag
    y = y .* hanning(length(y))';
end

wave = NaN(1,length(lfp_TimeStamps) * lfp_SamplesPerFrame);
if taperflag
    wave(1:lfp_time2index(lfp_Events(lfp_TrialIndex(1,1),1))) = ...
        lfp_Samples{filenum}(1:lfp_time2index( ...
        lfp_Events(lfp_TrialIndex(1,1),1) ));
end
trials = 1:size(lfp_TrialIndex, 1);
for trial = trials
    eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
    trialevents = lfp_Events(eventrange,:);
    reftime = trialevents( ...
        find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
        1 );
    if length(reftime) == 0
        warning('lfp_rmEP:noRef', ...
            'Skipping trial %d, no reference event', trial);
        continue
    else
        reftime = reftime(1);
    end
    startpoint = lfp_time2index(reftime + xstart);
    endpoint = min(startpoint + length(y) - 1, numel(wave));
    if any(~isnan(wave(startpoint:endpoint)))
        warning('lfp_rmEP:overlap', ...
            'The time interval for trial %d overlaps with a previous trial', ...
            trial );
    end
    wave(startpoint:endpoint) = ...
        lfp_Samples{filenum}(startpoint:endpoint) - ...
        y(1:endpoint-startpoint+1);
    if taperflag
        if trial < size(lfp_TrialIndex, 1)
            idx2 = lfp_time2index(lfp_Events(lfp_TrialIndex(trial+1,1),1));
        else
            idx2 = numel(lfp_Samples{filenum});
        end
        idx1 = lfp_time2index(lfp_Events(lfp_TrialIndex(trial,1),1));
        wave(idx1:(startpoint-1)) = ...
            lfp_Samples{filenum}(idx1:(startpoint-1));
        wave((endpoint+1):idx2) = ...
            lfp_Samples{filenum}((endpoint+1):idx2);
    end
end
units = sprintf('%s', lfp_SamplesUnits{filenum});

