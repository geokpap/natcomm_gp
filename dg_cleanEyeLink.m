function samples = dg_cleanEyeLink(samples, Fs)
%INPUTS
% samples: EyeLink time series.
% Fs: sampling rate in Hz.
%OUTPUTS
% samples: same as input, but with fast blink-like glitches removed.
%NOTES
% Assumes that <samples> is in volts using the same output scale as in the
% prototype Georgios session of 'Prez/ApAv/2023-03-03_11-14-11'.  Computes
% <blinkfrac> using a 90 ms moving window.  We designate any region where
% <blinkfrac> is greater than 0 as a contaminated region, unless the
% region actually touches 1, in which case we consider it to be a
% legitimate blink (provided that it isn't too close to a run of glitches).

%$Rev:  $
%$Date:  $
%$Author: dgibson $

% <windur> is for computing <blinkfrac>, the fraction of values that are at
% "blink level", i.e. below <blinklev>.
windur = 0.09; 
blinklev = -0.05;
delval = -0.059; % to substitute for deleted values

isblink = samples < blinklev;
winsamps = round(Fs * windur);
if ~mod(winsamps, 2)
    winsamps = winsamps + 1;
end
blinkfrac = conv(single(isblink(:)), ones(1, winsamps), 'same') / winsamps;
isblinky = blinkfrac ~= 0;
% <startidx>, <endidx> represent start and end indices of candidate
% deletion regions:
[startidx, endidx] = dg_findruns(reshape(isblinky, 1, []));
% If any deletion regions run off the end, forget about them: 
if endidx(end) == numel(samples)
    startidx(end) = [];
    endidx(end) = [];
end
if startidx(1) == 1
    startidx(1) = [];
    endidx(1) = [];
end
samp2del = false(numel(samples), 1);
for regidx = 1:length(startidx)
% Now we check for legitimate blinks, which we will not delete:
    isgoodblink = any(blinkfrac(startidx(regidx):endidx(regidx)) == 1);
    if isgoodblink
        % find "edge samples" around the first run of 1s, i.e. the sample
        % before the first 1 value, and the first non-1 after the first 1.
        % <edgidx1>, <edgidx2> are indices into <blinkfrac>.
        edgidx1 = find( ...
            blinkfrac(startidx(regidx):endidx(regidx)) == 1, 1 ) - 1 ...
            + startidx(regidx) - 1;
        edgidx2 = find( blinkfrac(edgidx1 + 1 : endidx(regidx)) ~= 1, ...
            1 ) + edgidx1;
        hascleanedges = blinkfrac(edgidx1 - winsamps) == 0 ...
            && blinkfrac(edgidx2 + winsamps) == 0;
    else
        hascleanedges = false;
    end
    if ~(isgoodblink && hascleanedges)
        % Mark for deletion:
        samp2del(startidx(regidx):endidx(regidx)) = true;
    end
end
samples(samp2del) = delval;
