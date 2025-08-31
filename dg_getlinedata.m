function [evts, waves, extraevts] = dg_getlinedata(hF, nodetailsflag)
% <evts> is a column vector of event marker handles.
% <waves> is an array of line handles where row 1 points to red curves
% (normally upper limit), 2: blue curves (normally lower limit), 3: black
% curves (average). It is assumed that wave data will never contain two
% points, and any line containing two points is an event marker. <waves> is
% a 3-row rectangular array, but may contain trailing zeros if there are
% unequal numbers of lines of each color.  <extraevts> is a column vector
% of event marker handles that failed the "nodetailsflag" condition, or if
% <nodetailsflag> is false, it is empty.

%$Rev: 42 $
%$Date: 2009-10-10 01:00:18 -0400 (Sat, 10 Oct 2009) $
%$Author: dgibson $

evts = [];
waves = zeros(3, 0);
extraevts = [];
hA = get(hF, 'CurrentAxes');
hL = findobj(hA, 'Type', 'line');
xlen = zeros(1, length(hL));
for k = 1:length(hL)
    xlen(1, k) = length(get(hL(k), 'XData'));
end
isevtmarker = xlen == 2;
isgoodevtmarker = isevtmarker;
if nodetailsflag
    for k = find(isevtmarker)
        if isequal(get(hL(k), 'ButtonDownFcn'), ...
                'fprintf(1,''No details available\n'')' )
            if k < length(hL) ...
                    && isequal( get(hL(k+1), 'Color'), ...
                    get(hL(k), 'Color') ) ...
                    && all( ...
                    get(hL(k+1), 'XData') - get(hL(k), 'XData') < 1e-3 )
                set(hL(k), 'ButtonDownFcn', get(hL(k+1), 'ButtonDownFcn'));
            else
                x = get(hL(k), 'XData');
                warning('dg_comparePasteups:nodetails', ...
                    'The event following the nodetails event at %.6f fails time and color criteria', ...
                    x(1) );
            end
        else
            isgoodevtmarker(k) = false;
            extraevts(end+1,1) = hL(k);
        end
    end
end
evts = reshape(hL(isgoodevtmarker), [], 1);
ridx = 1;
bidx = 1;
kidx = 1;
for k = find(~isevtmarker)
    % Wave data
    linecolor = get(hL(k), 'Color');
    if isequal(linecolor, [1 0 0])
        waves(1, ridx) = hL(k);
        ridx = ridx + 1;
    elseif isequal(linecolor, [0 0 1])
        waves(2, bidx) = hL(k);
        bidx = bidx + 1;
    elseif isequal(linecolor, [0 0 0])
        waves(3, kidx) = hL(k);
        kidx = kidx + 1;
    else
        warning('getlinedata:color', ...
            'Unknown wave data color %s', dg_thing2str(linecolor) );
    end
end
end

