function [thresh, cohdisp] = lfp_coherogram3(labelling,d1,d2,...
    moving_win,trialave,pad,colorbarflag,...
    absflag,arrowflag,minspikes,NW,K,rmdcflag,rotateflag,p)
%lfp_coherogram3(labelling,data1,data2,moving_win,trialave,pad,...
%     absflag,arrowflag,minspikes,NW,K,rmdcflag,rotateflag,p)
% Based on Partha Mitra's chronux functions (http://chronux.org/).
% <data1>, <data2> are wave data (either or both may be a spike train
% presented as wave data, i.e. 1 at the sample closest to each spike and
% zero elsewhere). <data1>, <data2> are in channels x time form, i.e. each
% wave is one row. OR, they are same-length cell arrays of data in that
% format, in which case the coherogram is computed for each pair of cells
% (one from data1 and one from data2) and then averaged.  The number of
% samples in each cell must be the same, although the number of channels
% may vary.  When averaging, treat NaN as 0.  If <absflag> is set,
% then absolute value is computed first before averaging, i.e. phase info
% is discarded.  If <arrowflag> is set, then phase arrows are
% overlaid on the coherogram, with magnitude = 1 at statistically
% significant points and magnitude = 0 elsewhere.  If <minspikes>
% is > 0, then <data2> is considered to be a spike train presented as wave
% data, and each coherogram window must have at least that number of spikes
% (totalled over all trials for each unit individually) in order to be
% displayed.  (Note that if <data2> contains any values other than 0 or 1,
% the assumptions of this 'minspikes' feature are violated.)  If there are
% 1 to <minspikes> - 1, then it is returned as Inf and shown as the maximal
% color value; if zero then it is returned as NaN and shown as the minimal
% color value. <NW>, <K> specify values of nw and k to Matlab's dpss
% function. If <rotateflag> is true, adds 90 degrees to the phase
% before drawing arrows, so zero phase is straight up instead of straight
% to the right.
%  When the arrows are pointing clockwise of zero (zero being horizontal to
%  the right, and clockwise being negative phase), that means that the
%  second channel is lagging the first.
% 
% NOTE: <trialave> ignored

%$Rev: 300 $
%$Date: 2013-05-01 18:47:38 -0400 (Wed, 01 May 2013) $
%$Author: dgibson $

lfp_declareGlobals;

% Ensure that data1 and data2 are matching cell arrays
if isequal(class(d1), 'cell')
    if ~isequal(class(d2), 'cell') || (numel(d2) ~= numel(d1))
        error('lfp_coherogram3:baddata', ...
            '<data1> and <data2> must be of same type and size' );
    end
    data1 = d1;
    data2 = d2;
else
    data1{1} = d1;
    data2{1} = d2;
end
for dataidx = 1:numel(data1)
    if (size(data1{dataidx}, 2) ~= size(data1{1}, 2)) || ...
            (size(data2{dataidx}, 2) ~= size(data1{1}, 2))
        error('lfp_coherogram3:baddata2', ...
            'All waves must contain same number of samples' );
    end
end

Ntr = 0;
Fs = 1/lfp_SamplePeriod;
if (NW < 1)
    % Comparing the values of NW = [1 1.5 2], the first dpss is closest to
    % Hamming for NW = 1.5 (or maybe 2 depending on your concerns).
    warning('lfp_coherogram3:lowNW', ...
        ['Time-bandwidth product = %d; values below 1 do not\n' ...
        'increase frequency resolution, but cause window shape\n' ...
        'to approach rectangular'], NW );
end
for dataidx = 1:numel(data1)
    tapers = [ NW K ];
    reftic=tic;
    [C,phi,t2,f] = ...
        dg_cohgramc(data1{dataidx}', data2{dataidx}', moving_win, ...
        tapers, pad, Fs, [0, Fs/2], ...
        0, trialave, rmdcflag );
    fprintf('In lfp_coherogram3\n');
    toc(reftic);
    fstart = 1;
    fend = length(f);
    if ~isempty(lfp_FreqLim)
        fstart = max(dg_binsearch(f, lfp_FreqLim(1)) - 1, 1);
        fend = min(dg_binsearch(f, lfp_FreqLim(2)), length(f));
        if f(fend - 1) == lfp_FreqLim(2)
            fend = fend - 1;
        end
        f = f(fstart:fend);
    end
    
    cohdisp(:,:,dataidx) = C(:,fstart:fend)';
    phidisp(:,:,dataidx) = phi(:,fstart:fend)';
    Ntr = Ntr + size(data1{dataidx},1);
end

% Some lines adapted from cohgramc
N1 = size(data1{1}, 2);
Nwin=round(Fs*moving_win(1)); % number of samples in window
Nstep=round(moving_win(2)*Fs); % number of samples to step through
winstart=[1:Nstep:N1-Nwin+1];

t = t2 + labelling.start;

if minspikes > 0
    for dataidx = 1:numel(data1)
        for win = 1:length(winstart)
            indx=winstart(win):winstart(win)+Nwin-1;    % from cohgramc
            nspikes = sum(sum(data2{dataidx}(:,indx)));
            if nspikes < minspikes && nspikes > 0
                cohdisp(:,win,dataidx) = Inf;
            end
        end
    end
end

if numel(data1) > 1
    cohdisp(find(isnan(cohdisp))) = 0;
    if absflag
        cohdisp = mean(cohdisp,3);
    else
        cohdisp = complex(cos(phidisp).*cohdisp, sin(phidisp).*cohdisp);
        cohdisp = abs(mean(cohdisp,3));
    end
    
end

myfields = fieldnames(labelling);
if ismember('hF', myfields)
    figure(labelling.hF);
end
hI = imagesc(t,f,cohdisp);
hA = get(hI,'Parent');
set(hA, 'YDir', 'normal');
if ismember('xax', myfields)
    set(get(hA,'XLabel'),'String',labelling.xax);
else
    set(get(hA,'XLabel'),'String','Seconds');
end
if ismember('yax', myfields)
    set(get(hA,'YLabel'),'String',labelling.yax);
else
    set(get(hA,'YLabel'),'String','Hz');
end
if colorbarflag
    hCB = colorbar;
    set(get(hCB,'YLabel'), 'String', 'Coherence Magnitude');
end
if absflag
    thresh = tinv(1-p,K*Ntr/2)/sqrt(K*Ntr/2);   % seat-of-pants
    %               combination of Bijan's rule of thumb and his guess that
    %               probably taking abs cuts df in half; needs verification
    caveat = ' maybe';
else
    % The following formula is based on the approximation that the
    % distribution of coherence is Gaussian (Jarvis & Mitra 2001 section
    % 5.1) centered at zero (it is actually more like 0.05, from Jarvis &
    % Mitra 2001 eqn 5.4, but we normally look at coherence >> .05 so we
    % call it 0) with variance 1/sqrt(nu) where nu is the number of
    % independent observations (Mitra and Pesaran 1999).
    thresh = tinv(1-p,K*Ntr)/sqrt(K*Ntr);   % from Bijan's rule of thumb
    caveat = '';
end
statmsg = sprintf('; p=%g level: %g%s', p, thresh, caveat);
if ismember('title', myfields)
    title(labelling.title, 'Interpreter', 'none');
else
    title(sprintf('Coh nw=%d k=%d pad=%d%s', NW, K, pad, ...
        statmsg ), 'Interpreter', 'none');
end

if arrowflag
    if rotateflag
        phidisp = phidisp + pi/2;
    end
    arrowgridX = 1:ceil(size(cohdisp,2)/32):size(cohdisp,2);
    arrowgridY = 1:ceil(size(cohdisp,1)/24):size(cohdisp,1);
    arrowX = cos(phidisp(arrowgridY, arrowgridX));
    arrowY = sin(phidisp(arrowgridY, arrowgridX));
    signifmask = cohdisp(arrowgridY, arrowgridX) > thresh ...
        & ~isinf(cohdisp(arrowgridY, arrowgridX));
    arrowX = (arrowX .* signifmask);
    arrowY = (arrowY .* signifmask);
    Xgrid2 = 1:length(arrowgridX);
    Ygrid2 = 1:length(arrowgridY);
    Xmin = 1 - .5 * length(Xgrid2)/arrowgridX(end);
    Ymin = 1 - .5 * length(Ygrid2)/arrowgridY(end);
    Xmax = (size(cohdisp,2)+.5) * length(arrowgridX)/arrowgridX(end);
    Ymax = (size(cohdisp,1)+.5) * length(arrowgridY)/arrowgridY(end);
    
    hA2 = axes('Position', get(hA, 'Position'), 'Parent', get(hA, 'Parent'));
    set(gcf,'CurrentAxes',hA2);
    hQ = quiver(Xgrid2, Ygrid2, arrowX, arrowY, 0.6);
    set(hA2, ...
        'XLim', [Xmin, Xmax], ...
        'YLim', [Ymin, Ymax], ...
        'XTickLabel', '', ...
        'YTickLabel', '' ...
        );
    set(hA2, 'TickLength', [0 0]);
    set(hA2, 'Color', 'none');
    set(hQ, 'Color', 'k');
    set(gcf,'CurrentAxes',hA);  % so that later things have an imagesc to work on
end

