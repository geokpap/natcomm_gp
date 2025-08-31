function [nfiles, fnames, rates] = dg_sortBySampleRate(infile, outdir)
%dg_sortBySampleRate(infile, outdir)
% Breaks up a Neuralynx CSC file containing mixed sample rates into several
% output files each of which has a constant sample rate.  Saves results in
% dg_nlx2Mat format.  Will also read dg_nlx2Mat format files.
%INPUTS
% infile: absolute or relative pathname to a Neuralynx CSC file or
%   dg_nlx2Mat format CSC file.
% outdir: absolute or relative pathname to the directory into which output
%   files are to be written.
%OUTPUTS
% Actual data output all goes to files.  Filenames are constructed
%   automatically by appending to the name of <infile> a string of the form
%   '_f<N>' where <N> is a number specifying the sample rate in Hertz.  If a
%   file contains more than one section at a given sample rate, than a
%   letter is appended starting with a and running through the alphabet.
%   It is an error to try to process a file that contains more than 26
%   non-contiguous sections at the same sampling rate.
% nfiles: number of output files.
% fnames: names (relative to <outdir>) of output files.
% rates: column vector of sampling rates for each output file, same length
%   as <fnames>.
%NOTES
% The criterion for detecting a change in sample rate is that they differ
%   by more than one part in one thousand (0.1%).
% It is assumed that timestamps are given in microseconds for frames of 512
%   samples when calculating sample rates in Hertz.
% It is assumed that any pair of frames whose timestamps differ by more
%   than 1 second is bracketing a pause in recording rather than a sampling
%   rate change.
% It is assumed that the last frame is at the same sampling rate as the
%   second-to-last.

%$Rev: 230 $
%$Date: 2015-10-17 22:48:50 -0400 (Sat, 17 Oct 2015) $
%$Author: dgibson $

%   what was that function I wrote to handle 'ncs' vs 'mat' extensions
%   automatically?

[~, name, ext] = fileparts(infile);
switch lower(ext)
    case {'.ncs' '.dat'}
        [TS, Samples, header] = dg_readCSC(infile);
        dg_Nlx2Mat_SamplesUnits = 'AD'; %#ok<*NASGU>
        for k = 1:length(header)
            if regexp(header{k}, '^\s*-ADBitVolts\s+')
                ADBitVoltstr = regexprep(header{k}, ...
                    '^\s*-ADBitVolts\s+', '');
                ADBitVolts = str2double(ADBitVoltstr);
                if isempty(ADBitVolts)
                    warning('dg_Nlx2Mat:badADBitVolts', ...
                        'Could not convert number from:\n%s', ...
                        header{k} );
                else
                    Samples = ADBitVolts ...
                        * Samples;
                    dg_Nlx2Mat_SamplesUnits = 'V';
                end
            end
        end
    case '.mat'
        x = load(infile);
        if ~(ismember('dg_Nlx2Mat_Timestamps', fieldnames(x))) ...
                || ~(ismember('dg_Nlx2Mat_Samples', fieldnames(x)))
            error('dg_sortBySampleRate:badmat', ...
                '%s is not a dg_nlx2Mat format file.', infile);
        end
        TS = x.dg_Nlx2Mat_SamplesUnits;
        Samples = x.dg_Nlx2Mat_Samples;
        if ismember('dg_Nlx2Mat_SamplesUnits', fieldnames(x))
            dg_Nlx2Mat_SamplesUnits = x.dg_Nlx2Mat_SamplesUnits;
        else
            dg_Nlx2Mat_SamplesUnits = 'arbs';
        end
end

% Two-pass design: first find the breaks between sections of constant
% sample rate by finding the pairs of frames whose lengths differ by more
% than 0.1%, skipping over frames of length > 1 s.  Then create output
% files for each section of constant-length frames.  We assume that any
% skipped frames were actually the same duration as the preceding one,
% followed by a pause in recording, so as long as the next frame is still
% (almost) the same duration, it's part of the same run of constant(ish)
% sampling rate.
% First pass:
framedurs = reshape(diff(TS), [], 1);
diffs2skip = find(framedurs > 1e6);
origframenums = setdiff(1:length(framedurs), diffs2skip);
framedurs(diffs2skip) = [];
% ratios(k) is the ratio of the k+1st frame duration to the kth.
ratios = framedurs(2:end) ./ framedurs(1:end-1);
% <lastframes> contains indices into framedurs, which has had frames2skip
% deleted, and never did contain an entry for the very last (unmeasurable)
% frame.  Likewise, <ratios> is one shorter than <framedurs>.
lastframes = find(ratios > 1.001 | ratios < 1/1.001);
firstframes = [1; lastframes+1];
firstorgframes = [ origframenums(1)
    origframenums(lastframes + 1) ];
lastframes(end+1, 1) = length(framedurs);
sections = [firstframes lastframes];
nfiles = size(firstorgframes, 1);
fnames = cell(nfiles, 1);
rates = zeros(nfiles, 1); % sample rates in Hz
repnum = zeros(nfiles, 1); % number of times rates has repeated
ratesWithRepeats = [];
for secidx = 1:nfiles
    rates(secidx) = round(512e6 / median(framedurs( sections(secidx,1) : ...
        sections(secidx,2) )));
    if any(rates(1:secidx-1) == rates(secidx))
        ratesWithRepeats = union(ratesWithRepeats, rates(secidx));
        repnum(secidx) = sum(rates == rates(secidx));
    else
        repnum(secidx) = 1;
    end
end
% Second pass:
lastorgframes = [ firstorgframes(2:end) - 1
    length(TS)];
origsections = [firstorgframes lastorgframes];
% idiot check:
if length(TS) ~= sum(diff(origsections, 1, 2) + ...
        ones(nfiles, 1));
    error('dg_sortBySampleRate:idiot', ...
        'Idiot detected.');
end
for secidx = 1:nfiles
    if ismember(rates(secidx), ratesWithRepeats)
        suffix = 'a' + repnum - 1;
    else
        suffix = '';
    end
    fnames{secidx} = sprintf('%s_f%d%s.mat', name, rates(secidx), suffix);
    indices = origsections(secidx,1) : origsections(secidx,2);
    dg_Nlx2Mat_Timestamps = TS(indices);
    dg_Nlx2Mat_Samples = Samples(:, indices);
    save(fullfile(outdir, fnames{secidx}), 'dg_Nlx2Mat_Timestamps', ...
        'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
end

