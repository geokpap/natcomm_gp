function dg_VT2NlxEVMouse(vtfilename, evfilename, strobebit)
%dg_VT2NlxEV(vtfilename, evfilename)
%   Reads a VT events file ("Yasuo format" or "EACQ" file) from
%   <vtfilename> and writes it in Neuralynx Events file format to
%   <evfilename>.  "Does the right thing" with <strobebit> before writing.
%   The customary value for <strobebit> is 2^15, but for Takashi data it
%   should be 2^14.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

[TimeStamps, TTLIDs, trialends] = lfp_readYasuoEvents(vtfilename);
events(:,1) = 100 * TimeStamps;
events(:,2) = TTLIDs - strobebit;
outcomemarkers = repmat([0 1; 0 91; 0 2], length(trialends), 1);
outcomemarkers(:,2) = outcomemarkers(:,2) - strobebit;
outcomeidx = 1;
for endidx = trialends
    outcomemarkers(outcomeidx:outcomeidx+2, 1) = ...
        events(endidx,1) + [1.004e6; 1.5e6; 2e6];
    outcomeidx = outcomeidx+3;
end
events = sortrows([events; outcomemarkers]);
dg_writeNlxEvents(evfilename, events(:,1)', events(:,2)');
