function timestamp = lfp_index2time(index)
%timestamp = index2time(index)
% Computes the absolute timestamp from the index into the sample arrays.
% Cannot be used until lfp_SamplesPerFrame, lfp_TimeStamps and
% lfp_SamplePeriod have legitimate values in the same time units as
% lfp_TimeStamps.  'index' can be any array, and 'timestamp' will be an
% array of the same size.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_TimeStamps lfp_SamplePeriod lfp_SamplesPerFrame;

frame = fix((index - 1)/lfp_SamplesPerFrame) + 1;
timestamp = reshape(lfp_TimeStamps(frame), size(frame)) + ...
    lfp_SamplePeriod * (index - 1 - (frame - 1) * lfp_SamplesPerFrame);