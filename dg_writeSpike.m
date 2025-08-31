function dg_writeSpike(filename, TS, Samples, Hdr, varargin)
%dg_writeSpike(filename, TS, Samples, Hdr)
% Only works for Windows.  Writes in Single Electrode format if <Samples>
% is 2-D or has just a single column; Stereotrode format if <Samples> has
% two columns; Tetrode format if <Samples> has four columns.  Raises a
% warning if the extension of <filename> does not match the size of
% <Samples> (Single Electrode, *.NSE; Stereotrode, *.NST; Tetrode, *.NTT).
%INPUTS
% TS:  Row vector of timestamps in clock ticks (units as recorded).
% Samples:  waveform samples in (samples, wires, triggers) format, or
%   optionally in (samples, triggers) format for .NSE data.
% Hdr:  Neuralynx format header (<16kB); may be omitted if no options are
%   specified, or empty.
%OPTIONS
% 'cellnums', cellnums - specify cluster numbers (or "cell numbers", in
%   Neuralynx parlance).  <cellnums> is a vector of length size(Samples,3).
% 'nlxjunk', nlxjunk - <nlxjunk> is a struct as returned from
%   'dg_readSpike', and it is used to populate all the other fields that
%   for our usual purposes are junk.  This option is only available in
%   Windows.
%NOTES
%   The dg_write* series of functions for creating Neuralynx format files
%   only works on Windows.  The dg_save* series can be used to save .mat
%   files in lfp_lib-compatible format.


%$Rev: 279 $
%$Date: 2021-10-15 17:53:12 -0400 (Fri, 15 Oct 2021) $
%$Author: dgibson $

cellnums = [];
nlxjunk = [];

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'cellnums'
            argnum = argnum + 1;
            cellnums = varargin{argnum};
        case 'nlxjunk'
            argnum = argnum + 1;
            nlxjunk = varargin{argnum};
        otherwise
            error('dg_writeSpike:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(cellnums)
    if isempty(nlxjunk)
        cellnums = zeros(size(TS));
    else
        cellnums = nlxjunk.CellNumbers;
    end
else
    cellnums = reshape(cellnums, 1, []);
end

% NOTE: Mat2NlxSE_411 seems to have some sort of problem with header lines
% that contain more than 126 characters, and replaces characters number 127
% and 128 with a line break (CR/LF, I assume).
if nargin < 4 || isempty(Hdr)
    Hdr = {'######## Neuralynx Data File Header '
        sprintf('## File Name: %s ', filename)
        sprintf('## Time Opened: %s ', datestr(now))
        '## written by dg_writeSpike '
        ' '
        };
end

if size(Samples,1) ~= 32
    error('dg_writeSpike:Samples', 'There must be 32 Samples per trigger');
end

if length(size(Samples)) < 3 || size(Samples,2) == 1
    numwires = 1;
elseif size(Samples,2) == 2
    numwires = 2;
elseif size(Samples,2) == 4
    numwires = 4;
else
    error('dg_writeSpike:Samples2', 'Bad format for Samples');
end
[p, n, ext] = fileparts(filename);
switch upper(ext)
    case '.NSE'
        if numwires ~= 1
            warning('dg_writeSpike:ext', ...
                'File extension is .NSE but Samples contains %d wires', ...
                numwires);
        end
    case '.NST'
        if numwires ~= 2
            warning('dg_writeSpike:ext', ...
                'File extension is .NST but Samples contains %d wires', ...
                numwires);
        end
    case '.NTT'
        if numwires ~= 4
            warning('dg_writeSpike:ext', ...
                'File extension is .NTT but Samples contains %d wires', ...
                numwires);
        end
end
if ispc
    if isempty(nlxjunk)
        switch numwires
            case 1
                Mat2NlxSE_411(filename, 0, 1, 1, length(TS), ...
                    [1 0 1 0 1 1], TS, ...
                    cellnums, Samples, Hdr);
            case 2
                Mat2NlxTS_411(filename, 0, 1, 1, length(TS), ...
                    [1 0 1 0 1 1], TS, ...
                    cellnums, Samples, Hdr);
            case 4
                Mat2NlxTT_411(filename, 0, 1, 1, length(TS), ...
                    [1 0 1 0 1 1], TS, ...
                    cellnums, Samples, Hdr);
        end
    else
        switch numwires
            case 1
                Mat2NlxSE_411(filename, 0, 1, 1, length(TS), ...
                    [1 1 1 1 1 1], TS, nlxjunk.ScNumbers, cellnums, ...
                    nlxjunk.Params, Samples, Hdr);
            case 2
                Mat2NlxTS_411(filename, 0, 1, 1, length(TS), ...
                    [1 1 1 1 1 1], TS, nlxjunk.ScNumbers, cellnums, ...
                    nlxjunk.Params, Samples, Hdr);
                warning('dg_writeSpike:untested1', ...
                    'Untested code here.');
            case 4
                Mat2NlxTT_411(filename, 0, 1, 1, length(TS), ...
                    [1 1 1 1 1 1], TS, nlxjunk.ScNumbers, cellnums, ...
                    nlxjunk.Params, Samples, Hdr);
                warning('dg_writeSpike:untested2', ...
                    'Untested code here.');
        end
    end
elseif isunix
    switch numwires
        case 1
            Mat2NlxSE_v2(filename, 0, 1, 1, length(TS), [1 0 1 0 1 1], TS, ...
                cellnums, Samples, Hdr);
        case 2
            Mat2NlxTS_v2(filename, 0, 1, 1, length(TS), [1 0 1 0 1 1], TS, ...
                cellnums, Samples, Hdr);
        case 4
            Mat2NlxTT_v2(filename, 0, 1, 1, length(TS), [1 0 1 0 1 1], TS, ...
                cellnums, Samples, Hdr);
    end
else
    error('dg_writeSpike:arch', ...
        'Unrecognized architecture.');
end
