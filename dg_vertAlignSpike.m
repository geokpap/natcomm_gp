function dg_vertAlignSpike(infile, outfile)
% Shifts each waveform up or down so that the average of the pre-trigger
% samples is zero.  The "pre-trigger" sample range is defined as the first
% sample through the sample that is 1/2 of the way to the 'AlignmentPt'
% parameter saved in the Neuralynx file header.
%INPUTS
% infile: absolute or relative path to Neuralynx *.nse spike file.
% outfile: absolute or relative path to output file to be created; may be
%   either *.mat ('mcc' format, see NOTES) or *.nse (Neuralynx
%   Single Electrode format).
%OUTPUTS
% All outputs go to the output file.
%NOTES
% Currently only works on single-wire *.nse files.
% In *.nse files from the Cheetah "6.4.0 Development" version, the
%   'AlignmentPt' is actually the peak time; it seems the algorithm for
%   triggering and aligning a snapshot is to find the last sample that is
%   strictly increasing after the threshold has been exceeded.
% 'mcc' format has timestamps in seconds in col. 1, cluster  number in col.
%   2, and samples in the remaining columns.  This is the same format saved
%   by lfp_save_noGUI.

[TS, Samples, Hdr] = dg_readSpike(infile);
Samples = squeeze(Samples); % convert to (samples, triggers) format
% Extract 'AlignmentPt', save in <align>:
for k = 1:length(Hdr)
    if regexp(Hdr{k}, '^\s*-AlignmentPt\s+')
        AlignmentPtstr = regexprep(Hdr{k}, ...
            '^\s*-AlignmentPt\s+', '');
        align = str2double(AlignmentPtstr);
    end
end

% Do the shifting:
PTrange = 1 : round(align/2);
Samples = Samples - repmat( mean(Samples(PTrange, :), 1), ...
    size(Samples, 1), 1 );

% Save transformed data:
[~, ~, ext] = fileparts(outfile);
switch ext
    case '.mat'
        lfp_save_spikes = zeros(size(Samples,2), size(Samples,1) + 2);
        lfp_save_spikes(:, 1) = TS * 1e-6;
        lfp_save_spikes(:, 3:end) = Samples'; %#ok<NASGU>
        save(outfile, 'lfp_save_spikes');
    case '.nse'
        dg_writeSpike(outfile, TS, Samples, Hdr);
    otherwise
        error('dg_vertAlignSpike:unkext', ...
            'Unknown output filename extension "%s".', ext);
end

