function matrix = dg_xcovMatrix(sessiondir, filenames, varargin)
%INPUTS
% sessiondir: the directory from which files will be read.
% filenames: a cell array of filenames to read from <sessiondir>.
%   Filenames with extensions .ncs or .dat will be read as Neuralynx files,
%   files with .mat extensions will be read as dg_Nlx2Mat-format .mat
%   files.
%OUTPUTS
% matrix: the matrix of zero-lag cross-covariances between each pair of
%   files, in the same order as <filenames>.  The matrix is symmetrical by
%   definition.
%OPTIONS
% 'verbose': produces behavior of such a nature as might be considered to
%   justify the selection of the keyword 'verbose' for the purpose of
%   invoking this option.
%NOTES
% Based on spot-checking a few 

%$Rev: 211 $
%$Date: 2015-03-12 21:08:04 -0400 (Thu, 12 Mar 2015) $
%$Author: dgibson $

verboseflag = false;

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
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_xcovMatrix:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if verboseflag
    fprintf('Running in directory %s\n', sessiondir);
end
matrix = NaN(length(filenames));
for ix1 = 1:length(filenames)
    [TS1, Samples1] = readfile(fullfile(sessiondir, filenames{ix1}));
    for ix2 = ix1+1 : length(filenames)
        [TS2, Samples2] = readfile(fullfile(sessiondir, filenames{ix2}));
        if isequal(TS1, TS2)
            iscommon = true(size(TS1));
            iscommon2 = iscommon;
        else
            iscommon = ismember(TS1, TS2);
            if sum(iscommon) < length(TS1)/2
                error('dg_xcovMatrix:badTS', ...
                    'There is only %d%% overlap between timestamps', ...
                    round(100 * sum(iscommon) / length(TS1)) );
            end
            iscommon2 = ismember(TS2, TS1);
        end
        matrix(ix1, ix2) = xcov( ...
            reshape(Samples1(:,iscommon), [], 1), ...
            reshape(Samples2(:,iscommon2), [], 1), 0, 'coeff' );
        matrix(ix2, ix1) = matrix(ix1, ix2);
        if verboseflag
            fprintf('Finished %s vs %s\n', filenames{ix1}, filenames{ix2});
        end
    end
end
end

function [ts, samp] = readfile(filename)
[~,~,ext] = fileparts(filename);
switch lower(ext)
    case {'.ncs' '.dat'}
        [ts, samp] = dg_readCSC(filename);
    case {'.mat'}
        s = load(filename);
        ts = s.dg_Nlx2Mat_Timestamps;
        samp = s.dg_Nlx2Mat_Samples;
    otherwise
        error('dg_xcovMatrix:ext', ...
            'Unknown file extension: %s', ext);
end
end

