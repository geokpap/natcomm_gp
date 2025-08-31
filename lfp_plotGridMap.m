function hF = lfp_plotGridMap(varargin)
%hF = lfp_plotGridMap(trials, filenums, window, coords)
% Given a vector of filenums and a 2-column array of (x,y) coordinates
% corresponding to each filenum, pseudocolor plot the value of each filenum
% averaged over <window>.  Default plotting style is a simple grid whose
% spacing is median(diff(unique(coords(:,1)))) in the horizontal direction
% and median(diff(unique(coords(:,2)))) in the vertical direction.
%INPUTS
% trials, filenums, window: lfp_lib standard, per lfp_CSCboilerplate,
%   except that <trials> should be a scalar.
% coords: 2-column array of (x,y) coordinates with one row for each element
%   of <filenums>, col. 1 = x, col. 2 = y.
%OUTPUTS
% hF: handle to newly plotted figure.
%OPTIONS
% All options that can be given to lfp_getSamples can be given to
% lfp_disp2d with the same effect.  In addition, the following options are
% available:
% 'axes', hA - plots into the pre-existing axes handle <hA>.
% 'clim', clim - instead of setting the color scale to span the range of
%   the data, the range from clim(1) to clim(2) is used.
% 'dB' - transforms CSC values to dB (10*log10(v)) before plotting.
% 'equal' - sets axes aspect ratio using axis('equal').
% 'movie', avifilename - create a Matlab movie and save it as a 25
%   frames/sec (same rate as PAL/SECAM) AVI file at the location specified
%   by <avifilename>, which can be an absolute or relative path to the
%   output file.  Instead of averaging over all of <window>, the average
%   for each frame is computed over each 40 ms sub-window that is
%   completely contained within <window>, starting at <window(1)>.
%   Heedlessly blows away any pre-existing file at that location.  Do not
%   touch the figure window that is creating the movie until it is
%   finished.
% 'movtimescale', scale - magnifies the time scale of a movie made using
%   the 'movie' option by a factor of <scale> (e.g. scale = 2 means that
%   the movie will take twice as long to play as it took to record the
%   data).  <scale> can be any positive floating-point number.
% 'norm' - normalizes all values relative to the highest value in the raw
%   (i.e. before averaging) data.
% 'xlabel', xlabelstr - x axis label string (default is none).
% 'ylabel', ylabelstr - y axis label string (default is none).

%$Rev: 379 $
%$Date: 2016-03-30 15:21:20 -0400 (Wed, 30 Mar 2016) $
%$Author: dgibson $

global lfp_SamplePeriod lfp_SamplesUnits

[trials, filenums, window, arglist, verboseflag, getSamplesOpts] = ...
    lfp_CSCboilerplate(varargin); 
[sampledata, ~, ~, ~, badtrials, trials, filenums] = ...
    lfp_getSamples(trials, filenums, window, getSamplesOpts{:});
trials = setdiff(trials, badtrials);
if isempty(getSamplesOpts) && isempty(arglist)
    optstring = '';
else
    optstring = dg_thing2str([getSamplesOpts arglist]);
end
coords = arglist{1};
if ~isscalar(trials)
    error('lfp_plotGridMap:trials', ...
        '<trials> must contain just one trialnum.');
end
if size(coords,1) ~= numel(filenums)
    error('lfp_plotGridMap:coords', ...
        '<coords> must contain one row for each filenum.');
end

arglist(1) = [];
argnum = 1;
avifilename = '';
clim = [];
dBflag = false;
equalflag = false;
fps = 25;
hA = [];
mov = [];
normflag = false;
Nsub = 1;
scale = 1;
xlabelstr = '';
ylabelstr = '';
while argnum <= length(arglist)
    if ischar(arglist{argnum})
        switch arglist{argnum}
            case 'axes'
                argnum = argnum + 1;
                hA = arglist{argnum};
            case 'clim'
                argnum = argnum + 1;
                clim = arglist{argnum};
            case 'dB'
                dBflag = true;
            case 'equal'
                equalflag = true;
            case 'movie'
                argnum = argnum + 1;
                avifilename = arglist{argnum};
            case 'movtimescale'
                argnum = argnum + 1;
                scale = arglist{argnum};
            case 'norm'
                normflag = true;
            case 'xlabel'
                argnum = argnum + 1;
                xlabelstr = arglist{argnum};
            case 'ylabel'
                argnum = argnum + 1;
                ylabelstr = arglist{argnum};
            otherwise
                error('lfp_plotGridMap:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(arglist{argnum}));
        end
    else
        error('lfp_plotGridMap:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(arglist{argnum}));
    end
    argnum = argnum + 1;
end

% Set grid:
xmin = min(coords(:,1));
xmax = max(coords(:,1));
xbin = median(diff(unique(coords(:,1))));
% nudge xmax a little higher so that the actual max value will be strictly
% less than xmax:
xmax = xmax + xbin / 256;
numxbins = ceil((xmax - xmin) / xbin);
xmargin = (numxbins * xbin - (xmax - xmin)) / 2;
x0 = xmin - xmargin;
xedges = x0 + (0:numxbins) * xbin;
ymin = min(coords(:,2));
ymax = max(coords(:,2));
ybin = median(diff(unique(coords(:,2))));
ymax = ymax + ybin / 256;
numybins = ceil((ymax - ymin) / ybin);
ymargin = (numybins * ybin - (ymax - ymin)) / 2;
y0 = ymin - ymargin;
yedges = y0 + (0:numybins) * ybin;

% convert <coords> to grid indices <grididx>:
grididx = zeros(size(coords));
for coordsrow = 1:size(coords,1)
    grididx(coordsrow, :) = [ find( xedges(1:end-1) <= coords(coordsrow,1) ...
        & xedges(2:end) > coords(coordsrow,1) ) ...
        find( coords(coordsrow,2) >= yedges(1:end-1) ...
        & coords(coordsrow,2) < yedges(2:end) ) ];
end

if ~isempty(avifilename)
    framesize = round(1 / (Nsub * scale * lfp_SamplePeriod * fps));
end

if isempty(hA)
    hF = figure;
    hA = axes('Parent', hF);
else
    hF = get(hA, 'Parent');
end

% <sampledata> is assumed to contain just one trigger, so we squeeze away
% the second dimension, leaving us with sampledata(samples, filenums):
sampledata = squeeze(sampledata);
if normflag
    sampledata = sampledata / max(sampledata(:));
end
if dBflag
    sampledata = 10 * log10(sampledata);
end
if isempty(clim)
    clim = [min(sampledata(:)) max(sampledata(:))];
end

if isempty(avifilename)
    hCB = plotoneframe(hA, sampledata, grididx, numybins, numxbins, ...
        clim, equalflag);
else
    % 'movie' has been invoked
    numframes = floor(size(sampledata, 1) / framesize);
    for framenum = 1:numframes
        range = ((framenum - 1) * framesize) + 1 : (framenum * framesize);
        hCB = plotoneframe( hA, sampledata(range, :), grididx, ...
            numybins, numxbins, clim, equalflag);
        if framenum == 1
            mov = getframe(hA);
        else
            mov(framenum) = getframe(hA); %#ok<AGROW>
        end
    end
    % Convert from Matlab "movie" to AVI file:
    if verboseflag
        fprintf('Creating movie file %s\n', avifilename);
    end
    if isempty(mov)
        error('lfp_plotGridMap:movie', ...
            'The movie is empty.');
    end
    v = ver('matlab');
    vtok = regexp(v.Version, '^(\d+\.\d+)', 'tokens');
    if str2double(vtok{1}) >= 7.12
        vObj = VideoWriter(avifilename);
        open(vObj);
        writeVideo(vObj, mov);
        close(vObj);
    else
        aviobj = avifile(avifilename, 'fps', fps); %#ok<DAVIFL>
        for framenum = 1:length(mov)
            aviobj = addframe(aviobj, mov(framenum));
        end
        aviobj = close(aviobj); %#ok<NASGU>
    end
end
if dBflag
    CBlabelstring = 'dB';
else
    CBlabelstring = lfp_SamplesUnits{filenums(1)};
end
labeleverything(hA, hCB, trials, filenums, window, optstring, ...
    xlabelstr, ylabelstr, CBlabelstring);
end


function hCB = plotoneframe(hA, sampledata, grididx, numybins, ...
    numxbins, clim, equalflag)
matrix = NaN(numybins, numxbins);
for fnidx = 1:size(sampledata,2)
    matrix(grididx(fnidx,2), grididx(fnidx,1)) = ...
        nanmean(sampledata(:, fnidx));
end
imagesc(matrix, 'Parent', hA);
caxis(hA, clim);
if equalflag
    axis(hA, 'equal');
    set(hA, 'Color', 'none');
end
hCB = colorbar('peer', hA);
end


function labeleverything(hA, hCB, trials, filenums, window, optstring, ...
    xlabelstr, ylabelstr, CBlabelstring)
global lfp_FileNames
ylabel(hCB, CBlabelstring);
if ~isempty(xlabelstr)
    xlabel(hA, xlabelstr);
end
if ~isempty(ylabelstr)
    ylabel(hA, ylabelstr);
end
clickstr = sprintf('Files:\n');
for fnidx = 1:length(filenums)
    clickstr = sprintf('%s%s\n', clickstr, lfp_FileNames{filenums(fnidx)});
end
clickstr = sprintf('%sOptions:\n%s', clickstr, optstring);
lfp_createFigTitle(hA, 'Grid Map', trials, window, ...
    'click for details', clickstr);
end



