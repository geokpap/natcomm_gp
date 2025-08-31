function lfp_plotFishAnalysisFile(filepath, varargin)
%INPUTS
% filepath: absolute or relative pathname to a file containing results
%   as saved by lfp_fishAnalysis, i.e. with values for the following
%   variables:
%       'EPdata', 'extrema', 'opts', 'files', 'isgoodfile', 'blockdirs'
%OPTIONS
% 'refidx', refidx - plots only the specified channel group(s).
% 'rowidx', rowidx - for each channel group plotted, plots only the
%   specified row(s), where rows are numbered sequentially across figures.

%$Rev: 289 $
%$Date: 2012-12-10 17:42:50 -0500 (Mon, 10 Dec 2012) $
%$Author: dgibson $

rows2plot = [];
refs2plot = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'refidx'
            argnum = argnum + 1;
            refs2plot = varargin{argnum};
        case 'rowidx'
            argnum = argnum + 1;
            rows2plot = varargin{argnum};
        otherwise
            error('lfp_plotFishAnalysisFile:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

load(filepath, '-mat');
if ~isempty(refs2plot)
    EPdata = EPdata(refs2plot);
    extrema = extrema(refs2plot);
    files = files(refs2plot);
end
for refidx = 1:length(files) %#ok<*USENS>
    filenames{refidx} = {files{refidx}.name}; %#ok<AGROW>
    if ~isempty(rows2plot)
        EPdata{refidx} = EPdata{refidx}(:,:,rows2plot,:);
        extrema{refidx} = extrema{refidx}(rows2plot,:);
        filenames{refidx} = filenames{refidx}(rows2plot);
    end
end
lfp_plotFishAnalysis(EPdata, extrema, filenames, opts, blockdirs);
