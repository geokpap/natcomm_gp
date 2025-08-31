function [badtrials, analysis] = lfp_findbadtrials(varargin)
%badtrials = lfp_findbadtrials(filenum)
%badtrials = lfp_findbadtrials(filename, sessiondir)
%   Applies a simple heuristic to identify trials in <filename> that are
% likely to be corrupted by some bizarre, enormous transient that should be
% excluded from analysis: computes the absolute difference between each
% sample and the mean for the whole channel, and returns the trialnums of
% any trials where that absolute difference exceeds 20 times the median
% value for the whole channel (i.e. any trial that contains a "bad value").
% Mean and median values are computed across all samples, unless 'windows'
% is specified.  For final analyses, you will want to specify 'windows'.
% The only advantage of NOT specifying it is running speed.  However, since
% lfp_time2trial considers the ITI following a given trial to belong to
% that trial, there may be trials listed in <badtrials> that only have bad
% values during the following ITI if 'windows' is not used.
%INPUTS
% filename|filenum: when called without specifying <sessiondir>, this
%   should be a filenum for a channel that is already loaded.  When
%   <sessiondir> is also given, it should be the name of a CSC file
%   (Neuralynx or .mat) in <sessiondir>.
% sessiondir: (not used if <filename> is a number) The path to a directory
%   from which to read session data using lfp_read2.
%OUTPUTS
% badtrials: trialnums, suitable for assignment to lfp_BadTrials or to use
%   as the first argument, <trials>, in a typical analysis function call.
% analysis: contains various results from the signal analyses that are used
%   to identify bad trials.  Struct with fields:
%       flatliners: trialnums that were marked bad because they contained
%         stretches of constant value; this is a subset of <badtrials>.
%       hardclip:  <true> if there is hard clipping, <false> if not.  Empty
%         if 'softclipthresh' is not used. This flag is only set if at
%         least 1 in 1000 samples is clipped.
%       medmag: median magnitude of selected samples after removing DC.
%       peakvals: sample magnitudes at which there were peaks in the sample
%         histogram.  Empty if 'softclipthresh' is not used.
%       softclipthresh: threshold determined from soft clipping
%         identification, or empty if there was no soft clipping
%         identified.
%       thresh: the actual threshold used to detect bad trials.
%OPTIONS
% 'flatliners', numpts - searches for sections where the same sample value
%   occurs more than <numpts> in a row, and marks any trial containing one
%   as bad. Has no effect if 'windows' is not specified.
% 'magfactor', magfactor - This value overrides the factor of 20 that is
%   used as a threshold for identifying bad trials.
% 'softclip' - checks for soft clipping, i.e. a range of values of
%   magnitude well below the actual hard clipping limits, but which
%   nonetheless produce a narrow peak (at least 2x lower sigma than the
%   main peak) in the sample histogram which is out in the tail of
%   distribution (beyond the 90th prctile).  This can easily happen in CSD
%   or other differential derivations.  Has no effect if 'windows' is not
%   specified.
% 'trial' - in addition to whatever 'windows' are specified (see 'windows'
%   option), also scans the entire trial as a whole.  Has no effect if
%   'windows' is not specified.
% 'windows', aligns, offsets - instead of searching the entire trial and
%   only the trial for bad values, searches a series of specified analysis
%   windows that are defined analogously to the analysis window specified
%   by the third argument for most of the lfp_lib functions.  <aligns> is a
%   cell vector of event IDs to use as alignment events, and <offsets> is a
%   two-column numeric array containing one row for each element of
%   <aligns> where the first column specifies the window start time
%   relative to the alignment event and the second column species the
%   window end time.  Each element in <aligns> is used together with its
%   row in <offsets> to determine a range of samples to search for bad
%   values.  Any trial that has bad values in any of the windows will be
%   returned in the <badtrials> list.  Also, when computing the median and
%   mean sample magnitudes, only data from enabled trials (see
%   lfp_enabledTrials) are used.
%EXAMPLES
% badtrials = lfp_findbadtrials(1);
%   Searches filenum 1 for trials containing bad values.
% badtrials = lfp_findbadtrials('lfp1.dat', '/smbshare/rodent/v1/s36/acq10')
%   Loads and then searches '/smbshare/rodent/v1/s36/acq10/lfp1.dat' for
%   trials containing bad values.
% badtrials = lfp_findbadtrials(1, 'windows', {2}, [5 7])
%   Search each trial that contains an event 2 for bad values that are
%   between 5 and 7 seconds after the event 2.
% badtrials = lfp_findbadtrials(1, 'windows', {2; [124 127]}, [-1 0; 1.3 3.4])
%   Search each trial that contains an event 2 for bad values that are
%   within the last second before event 2, and trials containing either of
%   events 124 or 127 for bad values between 1.3 and 3.4 s after the found
%   event (if events 124 and 127 are both found, the earlier one is used).
% badtrials = lfp_findbadtrials(1, 'windows', {[124 127]}, [0 3.5], 'trial');
%   Search the entirety of each trial for bad values, plus trials
%   containing event 124 or 127 for bad values within the first 3.5 s after
%   the event.  Note that this only makes sense to do if there are trials
%   that end earlier than 3.5 s after event 124 or 127.
%NOTES
%   See also lfp_getTrialID and lfp_getTrialNum for more information about
% "trialnums", "TrialIDs", and the relationship between them.

%$Rev: 404 $
%$Date: 2019-11-01 14:27:38 -0400 (Fri, 01 Nov 2019) $
%$Author: dgibson $

global lfp_Samples lfp_ActiveFilenums lfp_TrialIndex lfp_Events ...
    lfp_SelectedTrials lfp_BadTrials

filenum = [];
if isnumeric(varargin{1})
    filenum = varargin{1};
    varargin(1) = [];
else
    filename = varargin{1};
    sessiondir = varargin{2};
    varargin(1:2) = [];
end

aligns = {};
flatflag = false;
magfactor = 20;
numpts = 0;
softclipflag = false;
trialflag = false;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'flatliners'
            flatflag = true;
            argnum = argnum + 1;
            numpts = varargin{argnum};
        case 'magfactor'
            argnum = argnum + 1;
            magfactor = varargin{argnum};
        case 'softclip'
            softclipflag = true;
        case 'trial'
            trialflag = true;
        case 'windows'
            argnum = argnum + 1;
            aligns = varargin{argnum};
            argnum = argnum + 1;
            offsets = varargin{argnum};
        otherwise
            error('lfp_findbadtrials:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if flatflag && numpts < 2
    error('lfp_findbadtrials:numpts', ...
        '''flatliners'' <numpts> must be at least 2.');
end

if (flatflag || softclipflag) && isempty(aligns)
    warning('lfp_findbadtrials:options', ...
        'An option was specified that has no effect without ''windows''.');
end

if isempty(filenum)
    evtfname = dg_findEventsFile(sessiondir);
    lfp_read2('preset', sessiondir, {evtfname filename});
    filenum = lfp_ActiveFilenums(1);
end

analysis = [];

if isempty(aligns)
    mag = abs (lfp_Samples{filenum}(:) - mean(lfp_Samples{filenum}(:)) );
    analysis.medmag = median(mag);
    badidx = find(mag > magfactor * analysis.medmag);
    badtrials = unique(lfp_time2trial(lfp_index2time(badidx)));
    return
end

isenabledsample = true(size(lfp_Samples{filenum}));
disabledtrials = union(find(~lfp_SelectedTrials), lfp_BadTrials);
for trial = disabledtrials
    isenabledsample( ...
        lfp_TrialIndex(trial, 3) : lfp_TrialIndex(trial, 4) ) = false;
end
enabledsamples = lfp_Samples{filenum}(isenabledsample);
sampmean = mean(enabledsamples);

analysis.softclipthresh = []; % if empty, there was no soft clip peak found.
if softclipflag
    % Construct smoothed histogram of sample values and look for a side
    % peak whose width is less than half that of the main (i.e. at 0)
    % peak.
    absmag = abs(enabledsamples - sampmean);
    % make smoothed 1000-bin magnitude histo:
    numbins = 1000;
    binedges = linspace(0, max(absmag), numbins+1);
    binwidth = binedges(end) / numbins;
    binctrs = binedges + binwidth/2;
    counts = histc(absmag, binedges);
    smoothpts = 21;
    hw = hanning(smoothpts);
    hw = hw / sum(hw);
    smoocnts = conv(counts, hw, 'same');
    ispeak = [ false
        smoocnts(2:end-1) >= smoocnts(1:end-2) ...
        & smoocnts(3:end) < smoocnts(2:end-1)
        false ];
    analysis.peakvals = binctrs(ispeak);
    % Strictly speaking, there will be two distinct clipping values,
    % one positive and one negative.  However, how far apart they are
    % will depend on the value of '-ADBitVolts', which is not available
    % here.  We assume that they will form one smoothed peak that will
    % appear within <smoothpts> of the maximum magnitude:
    analysis.hardclip = analysis.peakvals(end) > binctrs(end-smoothpts) && ...
        sum(absmag > binctrs(end-(smoothpts))) ...
        / numel(absmag) >= 1e-3;
    % The whole rest of the analysis depends on there being at least
    % two peaks in the sample value histogram.  We designate the
    % biggest as the "main peak" and the second biggest as the "side
    % peak":
    if sum(ispeak) >= 2
        peakidx = find(ispeak);
        [~, peakidx2] = sort(smoocnts(peakidx), 'descend');
        % If the biggest peak is not the first peak, this whole model is
        % wrong, and the channel is probably bad:
        if peakidx2(1) ~= 1
            error('lfp_findbadtrials:mainpeak', ...
                'Bad channel - main peak is not the closest one to zero.');
        end
        sidepeakidx = peakidx(peakidx2(2));
        sidepeakval = binctrs(sidepeakidx);
        midpeakidx = round(mean([peakidx(1) sidepeakidx]));
        sidecnts = counts(sidepeakidx:end) - counts(midpeakidx);
        sidestd = binctrs(find( sidecnts < ...
            sidecnts(1) * normpdf(1)/normpdf(0), 1 ));
        analysis.softclipthresh = sidepeakval - 3 * sidestd;
    end
end
analysis.medmag = median(reshape( ...
    abs(lfp_Samples{filenum}(isenabledsample) - sampmean), [], 1 ));
% At this point we must adjudicate between <medmag> and <sidepeakval>:
if isempty(analysis.softclipthresh)
    % Nothing to decide!
    analysis.thresh = magfactor * analysis.medmag;
else
    % Choose the lower value:
    analysis.thresh = min(magfactor * analysis.medmag, analysis.softclipthresh);
end
isbadtrial = false(size(lfp_TrialIndex,1), 1);
isflatliner = false(size(lfp_TrialIndex,1), 1);
for trial = lfp_enabledTrials
    for alignidx = 1:numel(aligns)
        trialevts = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
        evtidx = find(ismember(trialevts(:,2), aligns{alignidx}), 1);
        if ~isempty(evtidx)
            evtTS = trialevts(evtidx, 1);
            winTS = evtTS + offsets(alignidx, :);
            winidx = lfp_time2index(winTS);
            trialmag = abs(lfp_Samples{filenum}(winidx(1):winidx(2)));
            if flatflag
                isconst = lfp_Samples{filenum}(winidx(1):(winidx(2)-1)) ...
                    == lfp_Samples{filenum}((winidx(1)+1):winidx(2));
                [startidx, endidx] = dg_findruns(isconst);
                runlen = endidx - startidx + 1;
                if any(runlen > numpts - 1)
                    isflatliner(trial) = true;
                end
            end
            if isflatliner(trial) || any(trialmag > analysis.thresh)
                isbadtrial(trial) = true;
                break
            end
        end
    end
    if ~isbadtrial(trial) && trialflag
        trialmag = abs(lfp_Samples{filenum}( ...
            lfp_TrialIndex(trial,3):lfp_TrialIndex(trial,4) ));
        if any( trialmag > analysis.thresh )
            isbadtrial(trial) = true;
        end
    end
end
badtrials = find(isbadtrial);
if flatflag
    analysis.flatliners = find(isflatliner);
end

