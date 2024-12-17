function sampleindex = lfp_time2index(time)
%sampleindex = lfp_time2index(time)
% Finds the sample index whose timestamp is closest to the specified
% <time>. Cannot be used until lfp_SamplesPerFrame and lfp_SamplePeriod
% have legitimate values in the same time units as lfp_TimeStamps. Although
% this works for time values that were not recorded, it could be a
% programming logic error to actually use it there.  <time> can be an
% array.  Return value is a guaranteed valid sample index.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_TimeStamps lfp_SamplePeriod lfp_SamplesPerFrame;

sampleindex = zeros(size(time));

for idx = 1:numel(time)
    nextframeindex = dg_binsearch(lfp_TimeStamps, time(idx));
    if nextframeindex > 1
        prevframeindex = nextframeindex - 1;
    else
        sampleindex(idx) = 1;
        continue
    end
    
    indexinframe = round( ...
        (time(idx) - lfp_TimeStamps(prevframeindex)) / lfp_SamplePeriod ) + 1;
    if indexinframe <= lfp_SamplesPerFrame
        sampleindex(idx) = indexinframe + ...
            lfp_SamplesPerFrame * (prevframeindex - 1);
        continue
    else
        % There was a gap between frames.
        endprevframe = lfp_TimeStamps(prevframeindex) + ...
            (lfp_SamplesPerFrame - 1) * lfp_SamplePeriod;
        if (nextframeindex > length(lfp_TimeStamps)) ...
                || ((time(idx) - endprevframe) < ...
                (lfp_TimeStamps(nextframeindex) - time(idx)))
            % return end of previous frame:
            sampleindex(idx) = lfp_SamplesPerFrame * prevframeindex;
        else
            % return beginning of next frame:
            sampleindex(idx) = lfp_SamplesPerFrame * prevframeindex + 1;
        end
    end
end
