function [result, units] = lfp_deglitch2(filenum, ...
    threshold, rt, falltime, endthresh, varargin)
%LFP_DEGLITCH2 function for use with lfp_createWave or lfp_createEvents.
%result = lfp_deglitch2(filenum, ...
%    threshold, risetime, falltime, endthresh)
%result = lfp_deglitch2(..., mode)
%result = lfp_deglitch2(..., 'markrejects')
%   Parametric deglitcher that offers more flexibility than lfp_deglitch.
%   Defines glitches as events that begin with a change in the sample
%   series <filenum> that exceeds <threshold> in magnitude, and whose sign
%   matches <threshold>, as compared to the value <risetime> samples
%   previous.  The start of the glitch is defined as the first sample
%   satisfying that condition, regardless of whether subsequent samples
%   also satisfy the condition.  The glitch is considered to be finished
%   when the signal is within <endthresh> of its pre-glitch value; however,
%   if it overshoots by more than the magnitude of <endthresh>, then the
%   entire glitch is skipped in confusion and reported with a warning and
%   log entry. The total duration of the glitch, defined as the number of
%   samples included between the threshold crossings, must be a minimum of
%   <falltime>(1) samples and a maximum of <falltime>(2) samples.  
%OPTIONS
% mode - If <mode> is 'evt', then a list of timestamps of the first sample
%   of each glitch is returned; if 'wave', a copy of the sample data in
%   <filenum> is returned with the glitch samples replaced by a linear
%   interpolation.  If 'both', then a cell array is returned containing
%   timestamps in result{1} and wave data in result{2}; note that in this
%   case, adapters are needed to interface lfp_deglitch2 with
%   lfp_createWave and lfp_createEvents. <mode> defaults to 'wave'.
% 'markrejects' - returns multiple event timestamp vectors in a cell
%   vector; first cell is standard list of glitch times, second is glitch
%   candidates that were rejected on the basis of the duration criterion,
%   third is glitch candidates that were rejected because they were
%   contained in another glitch, third is glitch candidate starts that end
%   with overshoot, fourth is glitch candidate ends that overshoot.
% 'bipolar' - sign of <threshold> does not matter, and glitches will be
%   marked in both directions.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

bipolar = false;
mode = 'wave';
markrejects = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'bipolar'
            bipolar = true;
        case {'both' 'evt' 'wave'}
            mode = varargin{argnum};
        case 'markrejects'
            markrejects = true;
        otherwise
            error('lfp_deglitch2:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isequal(size(filenum), [1 1])
    error('lfp_deglitch2:badfilenums', ...
        '<filenum> must be a single number' );
end
if ~isequal(size(threshold), [1 1])
    error('lfp_deglitch2:badthreshold', ...
        '<threshold> must be a single number' );
end
if ~isequal(size(rt), [1 1])
    error('lfp_deglitch2:badrisetime', ...
        '<risetime> must be a single number' );
end
if rt <= 0
    error('lfp_deglitch2:badrisetime2', ...
        '<risetime> must be greater than zero' );
end
if ~isequal(size(falltime), [1 2])
    error('lfp_deglitch2:badfalltime', ...
        '<falltime> must be a 1X2 array' );
end
if ~isequal(size(endthresh), [1 1])
    error('lfp_deglitch2:badendthresh', ...
        '<endthresh> must be a single number' );
end
if abs(endthresh) > abs(threshold)
    error('lfp_deglitch2:badendthresh2', ...
        '<endthresh> must not exceed <threshold> in magnitude' );
end
if endthresh ~= 0 && sign(endthresh) ~= sign(threshold)
    error('lfp_deglitch2:badendthresh3', ...
        '<endthresh> and <threshold> must have the same sign' );
end

numglitches = 0;
thresh = abs(threshold);
glitchlist = [];
wrongdur = [];
contained = [];
overshootstart = [];
overshootend = [];
markwrongdur = ismember(mode, {'evt', 'both'}) && markrejects;
if ismember(mode, {'wave' 'both'})
    wavedata = lfp_Samples{filenum};
end
units =lfp_SamplesUnits{filenum};

% Find diffs that cross threshold over rt samples.
trigpos = bipolar || threshold >= 0;
trigneg = bipolar || threshold < 0;
if trigneg
    negglitchidx = find(lfp_Samples{filenum}((1+rt):end) ...
        - lfp_Samples{filenum}(1:end-rt) < -thresh) + rt;
else
    negglitchidx = [];
end
if trigpos
    posglitchidx = find(lfp_Samples{filenum}((1+rt):end) ...
        - lfp_Samples{filenum}(1:end-rt) > thresh) + rt;
else
    posglitchidx = [];
end
glitchidx = union(posglitchidx, negglitchidx);

% For each glitchidx, iterate forwards over samples to find end of glitch.
% <refvalue> is the value at the last non-glitched sample.
thisentry = 1;
nsamples = numel(lfp_Samples{filenum});
while ~isempty(thisentry)
    startglitch = glitchidx(thisentry);
    trigwaspos = lfp_Samples{filenum}(startglitch-rt) < ...
        lfp_Samples{filenum}(startglitch);
    maxendglitch = falltime(2) + startglitch;
    endglitch = startglitch;
    refvalue = lfp_Samples{filenum}(startglitch-1);
    if trigpos
        posthresh = refvalue + endthresh;
    end
    if trigneg
        negthresh = refvalue - endthresh;
    end
    while (trigwaspos && ...
            lfp_Samples{filenum}(endglitch) > posthresh) ...
            || (~trigwaspos && ...
            lfp_Samples{filenum}(endglitch) < negthresh)
        endglitch = endglitch + 1;
        if endglitch > nsamples || endglitch > maxendglitch
            endglitch = Inf;
            break
        end
    end
    glitchdur = endglitch - startglitch;
    if glitchdur >= falltime(1) && glitchdur <= falltime(2)
        % Duration is good
        if abs(lfp_Samples{filenum}(endglitch) - refvalue) ...
                > abs(endthresh)
            % overshoot
            s = warning('query', 'lfp_deglitch2:overshoot');
            if isequal(s.state, 'on')
                msg = sprintf('%s glitch overshoot at sample %.0f, t=%.6f', ...
                    lfp_FileNames{filenum}, endglitch, lfp_index2time(endglitch) );
                warning('lfp_deglitch2:overshoot', '%s', msg );
                lfp_log(msg);
            end
            nextentry = thisentry + 1;
            if markrejects
                overshootstart(end+1) = startglitch;
                overshootend(end+1) = endglitch;
            end
        else
            % Glitch meets all criteria; process it
            numglitches = numglitches + 1;
            switch mode
                case 'evt'
                    glitchlist(end+1) = startglitch;
                case 'wave'
                    wavedata(startglitch-1 : endglitch) = ...
                        linspace(wavedata(startglitch-1), wavedata(endglitch), ...
                        endglitch - startglitch + 2 );
                case 'both'
                    glitchlist(end+1) = startglitch;
                    wavedata(startglitch-1 : endglitch) = ...
                        linspace(wavedata(startglitch-1), wavedata(endglitch), ...
                        endglitch - startglitch + 2 );
            end
            nextentry = find(glitchidx(thisentry+1:end) > endglitch);
            if ~isempty(nextentry)
                nextentry = nextentry(1) + thisentry;
            end
        end
    else
        % Duration no good
        if markwrongdur
            wrongdur(end+1) = startglitch;
        end
        msg = sprintf(...
            'Glitch at %.0f, t=%.6f, skipped because of duration %.0f', ...
            startglitch, lfp_index2time(startglitch), glitchdur );
        nextentry = thisentry + 1;
    end
    if nextentry > length(glitchidx)
        nextentry = [];
    end
    % Intervening entries in glitchidx are internal to the glitch; advance
    % to the next glitch start after the glitch end.
    if isempty(nextentry)
        if markrejects
            contained = [contained glitchidx(thisentry + 1 : end)];
        end
    else
        if markrejects
            contained = [contained ...
                glitchidx(thisentry + 1 : nextentry - 1) ];
        end
    end
    thisentry = nextentry;
end

if ismember(mode, {'evt', 'both'})
    if markrejects
        result{1} = {lfp_index2time(glitchlist)
            lfp_index2time(wrongdur)
            lfp_index2time(contained)
            lfp_index2time(overshootstart)
            lfp_index2time(overshootend)};
    else
        result{1} = lfp_index2time(glitchlist);
    end
end
switch mode
    case 'wave'
        result = wavedata;
    case 'both'
        result{2} = wavedata;
end
disp(sprintf('There were %d glitches.', numglitches ));
