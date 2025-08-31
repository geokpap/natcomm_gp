function varargout = lfp_bandcoh(trials, filenums, moving_win, ...
    freqlim, varargin)
%[hF, data] = lfp_bandcoh(trials, filenums, moving_win, freqlim)
%   Collects the common trial time in the usual manner (see
%   lfp_findCommonTime), computes coherence over trials and tapers in a
%   moving time window as in lfp_spec, computes the complex average
%   coherency over the frequency band f >= freqlim(1) & f <= freqlim(2),
%   and displays the average coherency magnitude and (optionally) phase as
%   in lfp_disp.
% INPUTS
%   trials - as usual
%   filenums - must be 2 element vector
%   moving_win - as for lfp_spec
%   freqlim - must be 2 element vector
% OUTPUTS
%   hF - handle to newly created figure
%   data - time points in first col, magnitude of coherency in second col,
%       and (optionally) phase in third column
% OPTIONS
%   Accepts 'nw', 'k', 'rmdc', 'rmEP', 'pad', as for lfp_mtspectrum.
%'err', [n p]
%  Uses chronux error curves for significance level <p>; <n>=1 for
%  theoretical calculation, <n>=2 for jackknife (takes longer).  WARNING:
%  the error intervals are calculated simply as the mean of the intervals
%  returned from Chronux across different frequencies.  I rationalize this
%  based on the approximation that the coherency values at all frequencies
%  in the <freqlim> band are 100% correlated.  IT IS THEREFORE AN
%  UNDERESTIMATE. However, using the sum of the squares of the error
%  intervals would produce a gross overestimate.  Someday (when we reach
%  the promised land) someone should reverse the order of the bootstrapping
%  and the averaging, i.e. rewrite the entire call stack from scratch. -DG
%'noplot'
%  Skips plotting and returns <data> as the first return value.
%'phi'
%  Displays phase in degrees as a second graph plotted below the first, as
%  in lfp_mtspectrum
%'phi-only'
%  Same as phi, but suppresses plotting of magnitude


%$Rev: 98 $
%$Date: 2009-11-22 21:29:26 -0500 (Sun, 22 Nov 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
if nargin < 3 || isempty(moving_win)
    moving_win = 1;
end
if size(moving_win,2) > 2 || ~(strcmp(class(moving_win), 'double'))
    error('lfp_bandcoh:badWindow2', ...
        '<moving_win> must contain one or two numbers');
end
if any(moving_win <= 0)
    error('lfp_bandcoh:badWindow3', ...
        '<moving_win> must contain positive, non-zero numbers');
end
if numel(moving_win) == 1
    moving_win(2) = moving_win/4;
end
if isempty(lfp_XLimAll)
    error('lfp_bandcoh:badxlim', ...
        'lfp_XLimAll must not be empty');
end

errflag = false;
K = 3;
NW = 2;
pad = 0;
phiflag = false;
magflag = true;
plotflag = true;
rmdcflag = false;
rmEPflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'err'
            errflag = true;
            argnum = argnum + 1;
            err = varargin{argnum};
        case 'k'
            argnum = argnum + 1;
            K = varargin{argnum};
        case 'noplot'
            plotflag = false;
        case 'nw'
            argnum = argnum + 1;
        case 'pad'
            argnum = argnum + 1;
            if ~strcmp(class(varargin{argnum}), 'double')
                error('lfp_bandcoh:badpad', ...
                    'The padding factor must be a number');
            end
            pad = varargin{argnum};
        case 'phi'
            phiflag = true;
        case 'phi-only'
            phiflag = true;
            magflag = false;
        case 'rmdc'
            rmdcflag = true;
        case 'rmEP'
            rmEPflag = true;
        otherwise
            error('lfp_bandcoh:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
tapers = [ NW K ];

% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

% Prepare data for chronux funcs in samples x trials form, using 3rd
% dimension for channels
oldxlimall = lfp_XLimAll;
lfp_XLimAll = [ lfp_XLimAll(1) - moving_win(1)/2 ...
    lfp_XLimAll(2) + moving_win(1)/2 ];
[interval, rawtrialinfo] = lfp_findCommonTime( ...
    lfp_enabledTrials(trials), 'recseg' );
if any(rawtrialinfo(:,3)==0)
    error('lfp_bandcoh:noref', ...
        'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(find(rawtrialinfo(:,3)==0))) );
end
xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
% Set the interval for analysis to xlimpoints clipped by <interval>:
interval(1) = max(xlimpoints(1), interval(1));
interval(2) = min(xlimpoints(2), interval(2));
idxrange = interval(1):interval(2);
% indices, data1, and data2 are all 2-D arrays of size [
% size(rawtrialinfo,1) x length(idxrange) ], i.e. trials x samples.
indices = ( repmat(idxrange, size(rawtrialinfo(:,3))) ...
    + repmat(rawtrialinfo(:,3), size(idxrange)) )';
if isempty(indices)
    error('lfp_mtspectrum:nodata', ...
        ['No samples were selected; note the if lfp_XLimAll is\n' ...
        'empty, <moving_win> is clipped to start and end of trial.'])
end
data1 = ...
    lfp_Samples{filenums(1)}(indices);
data2 = ...
    lfp_Samples{filenums(2)}(indices);
if rmdcflag
    data1 = data1 - repmat(mean(data1,1),size(data1,1),1);
    data2 = data2 - repmat(mean(data2,1),size(data2,1),1);
end
if rmEPflag
    data1 = data1 - repmat(mean(data1,2),1,size(data1,2));
    data2 = data2 - repmat(mean(data2,2),1,size(data2,2));
end
trialave = true;
Fs = 1/lfp_SamplePeriod;
% C and phi are time-windows x frequency index:
if errflag
    [C,phi,t,f,confC,phierr,Cerr] = ...
        dg_cohgramc(data1, data2, moving_win, ...
        tapers, pad, Fs, freqlim, err, trialave, rmdcflag );
else
    [C,phi,t,f] = ...
        dg_cohgramc(data1, data2, moving_win, ...
        tapers, pad, Fs, freqlim, 0, trialave, rmdcflag );
end
% t is in seconds re: the sample before data col. 1; data are taken from
% the sample range interval(1):interval(2)
t = t + (interval(1) + 1) * lfp_SamplePeriod;
% Converting from complex to polar and then back to complex and then back
% to polar is monstrously inefficient, but that commitment was made long
% ago in Cold Spring Harbor.  Note also that these conversions could be
% avoided if we wanted always to make the same assumption as is made by the
% code for the 'err' option, i.e. that the coherency is quasi-constant
% across <freqlim>.
S =  mean( (cos(phi) .* C) + 1i * (sin(phi) .* C), 2 );
data = [ t' abs(S) ];
if phiflag
    data(:, 3) = angle(S);
end

if plotflag
    hF = figure;
    varargout = {hF, data};
    if phiflag && magflag
        hA = subplot(2,1,1);
    else
        hA = axes;
    end
    if magflag
        plot(hA, t, data(:,2), 'k');
        set(hA, 'NextPlot', 'add');
        grid on;
        xlabel('Time, s');
        ylabel('Coherence');
        if errflag
            plot(hA, t, mean(Cerr(1,:,:),3)', 'r');
            plot(hA, t, mean(Cerr(2,:,:),3)', 'b');
        end
        if ~isempty(lfp_AlignmentRef)
            lfp_bandcoh_plotevtmarker(hA, lfp_AlignmentRef(1));
        end
    end
    if phiflag
        if magflag
            hA = subplot(2,1,2);
        end
        plot(hA, t, 180*data(:,3)/pi, 'k');
        set(hA, 'NextPlot', 'add');
        xlabel('Time, s');
        ylabel('Phase, deg');
        grid on;
        if errflag
            plot(hA, t, ...
                180*(data(:,3) + 2 * mean(phierr(1,:,:),3)')/pi, 'r');
            plot(hA, t, ...
                180*(data(:,3) - 2 * mean(phierr(1,:,:),3)')/pi, 'b');
        end
        if ~isempty(lfp_AlignmentRef)
            lfp_bandcoh_plotevtmarker(hA, lfp_AlignmentRef(1));
        end
    end
else
    varargout = {data};
end
lfp_XLimAll = oldxlimall;
end

function lfp_bandcoh_plotevtmarker(hA, evtid)
eventtime = 0;
lfp_declareGlobals;
eventcolor = '';
eventname = ''; % required for detailstr
if evtid <= length(lfp_EventNames)
    eventname = lfp_EventNames{evtid};
end
if evtid <= length(lfp_EventColors)
    eventcolor = lfp_EventColors{evtid};
else
end
if isempty(eventcolor)
    eventcolor = lfp_EventDefaultColor;
end
hL = plot([ eventtime eventtime ], ...
    get(hA, 'YLim'), ...
    'Color', eventcolor );
if isequal(class(eventcolor), 'char')
    eventcolorstr = eventcolor;
else
    eventcolorstr = mat2str(eventcolor);
end
detailstr = sprintf( ...
    '\\nAlignment Ref\\nLineColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
    eventcolorstr, eventname, evtid, evtid );
set(hL, ...
    'ButtonDownFcn', ...
    ['fprintf(1,''' detailstr '\n'')'] );
end
