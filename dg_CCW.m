function dg_CCW(infilename, outfilename, varargin)
%DG_CCW reads the Neuralynx tetrode data file <infilename>, applies the
%	cross-channel whitening of Emondi AA et al (2004) J Neurosci Methods
%	135:95, and writes the tetrode data with transformed waveforms to
%	<outfilename>.
%OPTIONS
%   'noCCW' - skips calculation of transform (to enable copying with gain
%       normalization).
%   'noclip' - deletes any spikes that have at least 3 samples in a row
%       with values 2047 or -2048 before gain normalization.
%   'normgain', mingain - normalizes gain of all 4 wires to mingain
%       (i.e. multiplies each channel by mingain/channelgain and
%       rounds).  If <mingain> is 0, then the gain of the lowest-gain
%       channel is used instead.  Note that this operation renders the
%       gains quoted in the Nlx header invalid in the output file.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

argnum = 1;
ccwflag = true;
noclipflag = false;
normflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'noCCW'
            ccwflag = false;
        case 'noclip'
            noclipflag = true;
        case 'normgain'
            normflag = true;
            argnum = argnum + 1;
            mingain = varargin{argnum};
    end
    argnum = argnum + 1;
end

[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints, NlxHeader] = ...
    Nlx2MatTT(infilename,1,1,1,1,1,1);

if noclipflag
    clipped = (DataPoints == 2047 | DataPoints == -2048);
    clipped3inrow = any(any( ...
        clipped(1:end-2,:,:) ...
        & clipped(2:end-1,:,:) ...
        & clipped(3:end,:,:), ...
        1), 2);
    TimeStamps(clipped3inrow) = [];
    ScNumbers(clipped3inrow) = [];
    CellNumbers(clipped3inrow) = [];
    Params(:, clipped3inrow) = [];
    DataPoints(:, :, clipped3inrow) = [];
end

if normflag
    ADgain = [];
    ampgain = [];
    for headnum = 1:length(NlxHeader)
        if regexp(NlxHeader{headnum}, '^\t-ADGain')
            ADgain = sscanf(NlxHeader{headnum}, ...
                '\t-ADGain\t%d\t%d\t%d\t%d');
        elseif regexp(NlxHeader{headnum}, '^\t-AmpGain')
            ampgain = sscanf(NlxHeader{headnum}, ...
                '\t-AmpGain\t%d\t%d\t%d\t%d');
        end
    end
    if isempty(ADgain) || isempty(ampgain)
        error('dg_CCW:missingheader', ...
            'ADgain=%s\nampgain=%s', ...
            dg_thing2str(ADgain), ...
            dg_thing2str(ampgain) );
    end
    channelgain = ampgain .* ADgain;
    if mingain == 0
        mingain = min(channelgain);
    end
    pointsperwire = size(DataPoints,2)/4;
    if fix(pointsperwire) ~= pointsperwire
        error('dg_CCW:badsamplesize', ...
            'Cannot divide samples in 4' );
    end
    DataPoints = round(DataPoints .* repmat(...
        mingain * ones(1, 4) ./ channelgain', ...
        [size(DataPoints,1) 1 size(DataPoints,3)] ));
end

if ccwflag
    % The "cross-channel correlation matrix" <CCcorr> is the correlation
    % between the entire set of samples from every spike on one wire <w1> with
    % the entire set of samples on another wire <w2>.  This matrix is
    % necessarily symmetric around the diagonal and has unity values along the
    % diagonal.
    for w = 1:4
        magnitudes(w) = sqrt(sum(sum(DataPoints(:, w, :) .^ 2)));
    end
    CCcorr = zeros(4);
    for w1 = 1:3
        for w2 = (w1 + 1) : 4
            CCcorr(w1, w2) = ...
                sum(sum(DataPoints(:, w1, :) .* DataPoints(:, w2, :), 1)) ...
                / (magnitudes(w1) * magnitudes(w2));
        end
    end
    % Fill in the missing values:
    CCcorr = CCcorr + CCcorr' + eye(length(CCcorr));

    % Find eigenvectors & eigenvalues:
    [V,D] = eig(CCcorr);
    lambda = diag(D,0);

    % Compute the transform matrix <Ecc>:
    Ecc = V' ./ repmat(sqrt(lambda), 1, 4);
else
    Ecc = eye(4);
end

% Transform each 4-channel sample to its cross-channel-whitened
% representation <OutputPoints>:
disp('Starting transform');
OutputPoints = zeros(size(DataPoints));
for m = 1:size(DataPoints,1)
    for n = 1:size(DataPoints,3)
        OutputPoints(m, :, n) = (Ecc * DataPoints(m, :, n)')';
    end
end

Mat2NlxTT(outfilename, 0, 1, 1, size(DataPoints,3), [1 1 1 1 1 1], ...
    TimeStamps, ScNumbers, CellNumbers, Params, OutputPoints, NlxHeader );

