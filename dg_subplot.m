function hA = dg_subplot(m, n, p, varargin)
% Much simpler than Matlab's subplot, which I decided not to use as
% starting material because it's too complicated and full-featured.
% hA = dg_subplot(m, n, p)
%   Creates a new axes at reading-order position <p> in <m> by <n> grid and
%   returns the handle to it.  All args are scalars.  As in "subplot", if
%   the axes already exists then nothing is created and the existing axes
%   becomes the current axes.
% hA = dg_subplot(m, n, [i j])
%   Creates a new axes at row <i>, column <j> in the grid.
%OPTIONS
% 'bottommargin', MarginSize - same as 'margin', but only applies to height
%   of margin at bottom of fig.  Defaults to ymargin.
% 'margin', MarginSize - specified in 'normalized' units.  Default value of
%   MarginSize is 0.1.
% 'xmargin', MarginSize - same as 'margin', but only applies to width of
%   margins.
% 'ymargin', MarginSize - same as 'margin', but only applies to height of
%   margins.
% 'paired' - no vertical space between plots below odd-numbered rows (or
%   equivalently, "above even-numbered rows")
% prop1, value1, prop2, value2, ... - as in "subplot", these are axes
%   properties and values for the new axes.
% 'topmargin', MarginSize - same as 'margin', but only applies to height
%   of margin at top of fig.  Defaults to ymargin.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

hF = gcf;

xmargin = 0.1;
ymargin = 0.1;
bottommargin = [];
topmargin = [];
pairedflag = false;
units = '';
argnum = 1;
args2delete = [];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'bottommargin'
            argnum = argnum + 1;
            bottommargin = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
        case 'margin'
            argnum = argnum + 1;
            xmargin = varargin{argnum};
            ymargin = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
        case 'paired'
            pairedflag = true;
            args2delete = [args2delete argnum];
        case 'topmargin'
            argnum = argnum + 1;
            topmargin = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
        case 'Units'
            argnum = argnum + 1;
            units = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
        case 'xmargin'
            argnum = argnum + 1;
            xmargin = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
        case 'ymargin'
            argnum = argnum + 1;
            ymargin = varargin{argnum};
            args2delete = [args2delete argnum-1 argnum];
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];
if isempty(bottommargin)
    bottommargin = ymargin;
end
if isempty(topmargin)
    topmargin = ymargin;
end


% Find I and J, respectively the row number and column number of the axes
% in the grid
if numel(p) == 1
    [J, I] = ind2sub([n m], p);
elseif numel(p) == 2
    I = p(1);
    J = p(2);
else
    error('dg_subplot:p', '<p> must be a scalar or two-element vector');
end

% Maintain grid as in "subplot" (more or less)
appdata = getappdata(hF);
appdatanames = fieldnames(appdata);
if ismember('SubplotGrid', appdatanames)
    grid = appdata.SubplotGrid;
else
    grid = NaN(m, n);
end

if isnan(grid(I, J))
    % Compute position in normalized units, which are used until the very
    % end, at which time they are converted to the value of the 'Units'
    % property (if given).
    width = (1 - (n + 1) * xmargin) / n;
    xpos = J * xmargin + max(0, (J - 1)) * width;
    if pairedflag
        if mod(m,2) == 0
            % even number of rows, use bottom margin
            height = ( 1 - ...
                (fix(m/2) - 1) * ymargin - topmargin - bottommargin ...
                ) / m;
            ypos = fix((m - I)/2) * ymargin + bottommargin ...
                + (m - I) * height;
        else
            % odd number of rows, no bottom margin
            height = (1 - (fix(m/2) - 1) * ymargin - topmargin)/m;
            ypos = fix((m - I + 1)/2) * ymargin + (m - I) * height;
        end
    else
        height = (1 - (m - 1) * ymargin - topmargin - bottommargin)/m;
        ypos = (m - I) * ymargin + bottommargin + (m - I) * height;
    end
    
    hA = axes('Position', [xpos ypos width height], varargin{:});
    if ~isempty(units)
        set(hA, 'Units', units);
    end
    grid(I, J) = hA;
    setappdata(hF, 'SubplotGrid', grid);
else
    % axes already exists
    hA = grid(I, J);
    set(hF, 'CurrentAxes', hA);
end
    
    

