function result = lfp_rmlinesc(filenum)
% DEPRECATED.  See lfp_rmlnspec.
%result = lfp_rmlinesc(filenum) for use with lfp_createWave.
% As far as I can tell, I abandoned this in favor of lfp_rmlnspec which
% both does more (removes harmonic series) and less (does not attempt to
% automatedly identify the frequencies to remove).  DG 11-Mar-2016.
% Also see dg_ch_rmlinesc (another deprecated function).

%  Wrapper for 8/24/2004 chronux rmlinesc.  Unfortunately, this uses a lot
%  of memory, so it must be done in batches.  Also, it is more effective
%  closer to the beginning of the waveform, so the actual processed
%  waveform extends beyond the end of each batch (if possible).  Frame
%  numbers of beginnings of batches are logged.  If possible, each trial is
%  processed as one batch so that there will not be splices in the middles
%  of trials.  Worse than all that, does not work very well.

%$Rev: 390 $
%$Date: 2017-06-21 17:37:26 -0400 (Wed, 21 Jun 2017) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(filenum), [1 1])
    error('lfp_rmlinesc:badfilenums', ...
        '<filenum> must have exactly 1 element' );
end

result = zeros(size(lfp_Samples{filenum}));
batchsize = 4;   % max number of frames
trialptr = 1;
frameptr = 1;
numframes = ceil(numel(lfp_Samples{filenum}/lfp_SamplesPerFrame));

while (numframes - frameptr + 1) > batchsize
    % There remains another batch of frames to process.
    % If possible, process through the next end of trial:
    sampleptr = (frameptr-1) * lfp_SamplesPerFrame + 1;
    currenttrial = find(lfp_TrialIndex(:,4) > sampleptr);
    if isempty(currenttrial)
        endsample = numel(lfp_Samples{filenum});
    else
        endsample = lfp_TrialIndex(currenttrial(1),4);
    end
    numframes2process = ceil( ...
        (endsample - sampleptr + 1) ...
        / lfp_SamplesPerFrame );
    if numframes2process <= batchsize
        lastframe = frameptr + numframes2process - 1;
    else
        lastframe = frameptr + batchsize - 1;
    end
    
    range = (sampleptr) : ...
        ((lastframe)*lfp_SamplesPerFrame);
    extension = round((lastframe - frameptr + 1)/2);
    extrange = (sampleptr) : ...
        min( numel(lfp_Samples{filenum}), ...
        ((lastframe + extension)*lfp_SamplesPerFrame) );
    if ~isempty(range)
        lfp_log(sprintf(...
            'lfp_rmlinesc starting batch starting with frame %d', ...
            frameptr));
        extresult = ...
            dg_ch_rmlinesc(lfp_Samples{filenum}(extrange)', ...
            [2 3], 1/lfp_SamplePeriod, [0 1/(2*lfp_SamplePeriod)], 3, .00005);
        result(range) = extresult(range - range(1) + 1);
        clear extresult;
        lfp_log(sprintf(...
            'lfp_rmlinesc finished batch starting with frame %d', ...
            frameptr));
    else
        break
    end
    frameptr = lastframe + 1;
end
% Now finish up the last fraction of a batch.
if ~isempty(range)
    range = (sampleptr) : ...
        numel(lfp_Samples{filenum});
    result(range) = ...
        dg_ch_rmlinesc(lfp_Samples{filenum}(range)', ...
        [2 3], 1/lfp_SamplePeriod, [0 1/(2*lfp_SamplePeriod)], 3);
end

