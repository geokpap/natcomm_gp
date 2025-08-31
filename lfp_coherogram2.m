function [thresh, p, cohdisp] = lfp_coherogram2(labelling,d1,d2,moving_win,trialave,pad,...
    absflag,arrowflag,minspikes)
%lfp_coherogram2(labelling,data1,data2,moving_win,trialave,pad,...
%     absflag,arrowflag,minspikes)
% <data1> is wave data, <data2> is a spike train presented as wave data,
% i.e. 1 at the sample closest to each spike and zero elsewhere.
% Plots based on Bijan Pesaran's functions <data1>, <data2> are in channels
% x time form, i.e. each wave is one row. OR, they are same-length cell
% arrays of data in that format, in which case the coherogram is computed
% for each pair of cells (one from data1 and one from data2) and then
% averaged.  The number of samples in each cell must be the same, although
% the number of channels may vary.  When averaging, treat NaN as 0.  If
% optional <absflag> is set, then absolute value is computed first before
% averaging, i.e. phase info is discarded.  If optional <arrowflag> is set,
% then phase arrows are overlaid on the coherogram, with magnitude = 1 at
% statistically significant points and magnitude = 0 elsewhere.  If
% optional <minspikes> is > 0, then each coherogram window must have at
% least that number of spikes (totalled over all trials for each unit
% individually) in order to be displayed.  If there are 1 to <minspikes> -
% 1, then it is returned as Inf and shown as the maximal color value; if
% zero then it is returned as NaN and shown as the minimal color value.
% 
% NOTE: <trialave> ignored, <pad> not implemented

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 9
    minspikes = 0;
end
if nargin < 8
    arrowflag = false;
end
if nargin < 7
    absflag = false;
end

% [coh_ptx, tfsp_X, tfsp_dN, rate, f] = ...
%     tfcoh_ptx(data1, data2, [moving_win(1) 5], 1/lfp_SamplePeriod, ...
%     moving_win(2), [0 1/lfp_SamplePeriod/2], pad, 0.05, trialave);
NW = 5; % desired time-bandwidth product
W = NW/moving_win(1);    % desired bandwidth in Hz

% Ensure that data1 and data2 are matching cell arrays
if isequal(class(d1), 'cell')
    if ~isequal(class(d2), 'cell') || (numel(d2) ~= numel(d1))
        error('lfp_coherogram2:baddata', ...
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
        error('lfp_coherogram2:baddata2', ...
            'All waves must contain same number of samples' );
    end
end

Ntr = 0;
for dataidx = 1:numel(data1)
    if trialave
        flag = 11;
    else
        flag = 11;
%         flag = 0;
    end
    [coh_ptx, tfsp_X, tfsp_dN, rate, f] = ...
        tfcoh_ptx(data1{dataidx}, data2{dataidx}, ...
        [moving_win(1) W], 1/lfp_SamplePeriod, moving_win(2), ...
        [0, 0.5/lfp_SamplePeriod], pad, 0.05, flag);
    
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
    
    cohdisp(:,:,dataidx) = coh_ptx(:,fstart:fend)';
    Ntr = Ntr + size(data1{dataidx},1);
end

% The range of samples in each window analyzed by tfcoh_ptx is given by
% X(:,win*dn+1:win*dn+n), where win is the window number starting with 1.
% So, the first window is centered at dn + (1+n)/2.
n = floor(moving_win(1)/lfp_SamplePeriod);
dn = floor(moving_win(2)/lfp_SamplePeriod);
nwin = size(coh_ptx, 1);
t = ((1+n)/2 + (1:nwin) * dn) * lfp_SamplePeriod;
t = t + labelling.start;

if minspikes > 0
    for dataidx = 1:numel(data1)
        for win = 1:nwin
            nspikes = sum(sum(data2{dataidx}(:,win*dn+1:win*dn+n)));
            if nspikes < minspikes && nspikes > 0
                cohdisp(:,win,dataidx) = Inf;
            end
        end
    end
end

if numel(data1) > 1
    cohdisp(find(isnan(cohdisp))) = 0;
    if absflag
        cohdisp = abs(cohdisp);
    end
    cohdisp = mean(cohdisp,3);
end

hI = imagesc(t,f,abs(cohdisp));
hA = get(hI,'Parent');
set(hA, 'YDir', 'normal');
colorbar;
K = floor(2*NW-1);  % number of tapers (see coherency_ptx)
p = .01;
if absflag
    thresh = tinv(1-p,K*Ntr/2)/sqrt(K*Ntr/2);   % seat-of-pants
    %               combination of Bijan's rule of thumb and his guess that
    %               probably taking abs cuts df in half
    caveat = ' maybe';
else
    thresh = tinv(1-p,K*Ntr)/sqrt(K*Ntr);   % from Bijan's rule of thumb
    caveat = '';
end
statmsg = sprintf('; p=%g level: %g%s', p, thresh, caveat);
title(sprintf('%s Coherogram%s', labelling.title, statmsg));

if arrowflag
    arrowgridX = 1:ceil(size(cohdisp,2)/50):size(cohdisp,2);
    arrowgridY = 1:ceil(size(cohdisp,1)/30):size(cohdisp,1);
    arrowX = real(cohdisp(arrowgridY, arrowgridX)) ...
        ./abs(cohdisp(arrowgridY, arrowgridX));
    arrowY = imag(cohdisp(arrowgridY, arrowgridX)) ...
        ./abs(cohdisp(arrowgridY, arrowgridX));
    signifmask = abs(cohdisp(arrowgridY, arrowgridX)) > thresh ...
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
    hQ = quiver(Xgrid2, Ygrid2, arrowX, arrowY, 0.8);
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

