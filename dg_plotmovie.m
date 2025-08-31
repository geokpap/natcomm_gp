function mov = dg_plotmovie(framesize, X, Y, varargin)
%mov = dg_plotmovie(framesize, X, Y, ...)
% End result is equivalent to
%   figure; plot(X,Y,...)
% but returns a Matlab movie of the process by recording a frame from
% the current axes every <framesize> points.  If either of X or Y is a
% matrix, then successive points are assumed to go down the columns of the
% matrix.  If the number of data points is not an integral multiple of
% <framesize>, the remaining points are plotted as the last frame.
%INPUTS
% framesize: number of data points per movie frame; need not be an integer.
% X, Y: data to submit to Matlab 'plot' function.
%OPTIONS
% All options and option values are passed through to the 'plot' function
% or the 'scatter' function verbatim, except for the following:
% 'append', mov - appends additional frames to the pre-existing movie
%   <mov>.
% 'axes', hA - plots into the pre-existing axes handle <hA>.
% 'scatter', S, C - uses the 'scatter' function in place of the 'plot'
%   function to take advantage of its marker coding features.  <S> is the
%   marker area argument and <C> is the marker color argument (see Matlab
%   help for 'scatter').  <X> and <Y> must be vectors when using 'scatter';
%   <S> and <C> can be scalars or vectors.
%NOTES
% It is extremely important not to do anything (e.g. move, resize, etc) to
% the figure window that is being used to record the movie because the
% movie recording mechanism and the graphics system easily get out of sync
% and cause frames to be captured that contain images of other windows.

%$Rev: 208 $
%$Date: 2014-10-29 17:41:36 -0400 (Wed, 29 Oct 2014) $
%$Author: dgibson $

if fix(framesize) ~= framesize
    warning('dg_plotmovie:framesize', ...
        'Rounding framesize to integer' );
    framesize = round(framesize);
end

if numel(X) == length(X)
    X = reshape(X, [], 1);
end
if numel(Y) == length(Y)
    Y = reshape(Y, [], 1);
end
if (numel(X) ~= length(X)) || (numel(Y) ~= length(Y))
    if size(X, 1) ~= size(Y, 1)
        error('dg_plotmovie:points', ...
            'The number of points (rows) in X and Y must be the same.' );
    end
end

argnum = 0;
C = [];
S = [];
hA = [];
offset = 0;
opts2delete = [];
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'append'
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            mov = varargin{argnum};
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            offset = length(mov);
        case 'axes'
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            hA = varargin{argnum};
            opts2delete(end+1) = argnum; %#ok<*AGROW>
        case 'scatter'
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            S = reshape(varargin{argnum}, [], 1);
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            C = reshape(varargin{argnum}, [], 1);
            opts2delete(end+1) = argnum; %#ok<*AGROW>
    end
end
varargin(opts2delete) = [];
if ischar(C)
    switch C
        case {'y' 'yellow'}
            C = [1 1 0];
        case {'m' 'magenta'}
            C = [1 0 1];
        case {'c' 'cyan'}
            C = [0 1 1];
        case {'r' 'red'}
            C = [1 0 0];
        case {'g' 'green'}
            C = [0 1 0];
        case {'b' 'blue'}
            C = [0 0 1];
        case {'w' 'white'}
            C = [1 1 1];
        case {'k' 'black'}
            C = [0 0 0];
        otherwise
            error('dg_plotmovie:color', ...
                'Unrecognized color string "%s".', C);
    end
end
if isequal(size(S), [1 1])
    S = repmat(S, size(X));
end
if isequal(size(C, 1), 1)
    C = repmat(C, size(X));
end

if isempty(hA)
    hF = figure;
    hA = axes('Parent', hF);
end
set(hA, 'NextPlot', 'add');

numframes = ceil(size(X, 1) / framesize);
mov(numframes + offset) = getframe(hA);
for framenum = 1:numframes
    lastpt = min(round(framenum * framesize), numel(X));
    if isempty(C)
        % standard line plot, requires joining to previously plotted
        % segment, which means we must back up by one point at the
        % beginning.
        firstpt = max(1, round((framenum-1) * framesize));
        plot(hA, X(firstpt:lastpt,:), Y(firstpt:lastpt,:), ...
            varargin{:});
    else
        firstpt = 1 + round((framenum-1) * framesize);
        scatter(hA, X(firstpt:lastpt,:), Y(firstpt:lastpt,:), ...
            S(firstpt:lastpt,:), C(firstpt:lastpt,:), varargin{:});
    end
    mov(framenum + offset) = getframe(hA);
end

