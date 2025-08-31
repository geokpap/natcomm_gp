function [result, firstoffender, lastoffender] = ...
    lfp_zeroxingFailure(filenum, startsample, endsample)
%OUTPUTS
% result: true if there are is a maximum with a negative value or a minimum
%     with a positive value in the interval from <startsample> to <endsample>
%     of the HHT waveform contained in <filenum>.  This is logically equivalent
%     (I think but have not formally proven) to failing the IMF criterion that
%     the number of zero crossings and the number of extrema differ by no more
%     than one.
% firstoffender: relative index into startsample:endsample of the first
%   zeroxing failure.
% lastoffender: relative index into startsample:endsample of the last
%   zeroxing failure.

%$Rev: 209 $
%$Date: 2011-02-19 16:13:00 -0500 (Sat, 19 Feb 2011) $
%$Author: dgibson $

global lfp_Samples;

firstoffender = [];
lastoffender = [];

relmaxima = dg_findpks(lfp_Samples{filenum}(startsample:endsample));
relminima = dg_findpks(-lfp_Samples{filenum}(startsample:endsample));

result = any(lfp_Samples{filenum}(relmaxima + startsample - 1) <= 0)  || ...
    any(lfp_Samples{filenum}(relminima + startsample - 1) >= 0);

if result
    offendingmins = ...
        find(lfp_Samples{filenum}(relminima + startsample - 1) >= 0);
    offendingmaxs = ...
        find(lfp_Samples{filenum}(relmaxima + startsample - 1) <= 0);
    if isempty(offendingmins)
        firstoffender = relmaxima(offendingmaxs(1));
        lastoffender = relmaxima(offendingmaxs(end));
    elseif isempty(offendingmaxs)
        firstoffender = relminima(offendingmins(1));
        lastoffender = relminima(offendingmins(end));
    else
        firstoffender = ...
            min([relminima(offendingmins(1)) relmaxima(offendingmaxs(1))]);
        lastoffender = ...
            max([relminima(offendingmins(end)) relmaxima(offendingmaxs(end))]);
    end
end