function barfolemew(data, GPcounts)
% Computes figures for "population.ai".
% GPcounts: numbers of units in classes [excitatory, inhibitory, 
%   unresponsive] according to GP (see Table p. 701,
%   "GeorgiosAnalysisVol5.docx"; the values in <GPcounts> can be copied
%   verbatim from GP email of 5/5/25, 9:40 AM).
numunresp = GPcounts(3);
numexcit = GPcounts(1);
firstE = numunresp + 1;
lastE = firstE + numexcit - 1;
baselines = mean(data(firstE:lastE, (-12:-1) + 41), 2);
for k = 1:numexcit
    renormed(k,:) = data(firstE + k - 1, :) - baselines(k); %#ok<AGROW>
end
p = (0.05/38)^(1/3);
thresh = -norminv(p);
isbig = abs(renormed) > thresh;
is3inrow = isbig(:, 1:end-2) & isbig(:, 2:end-1) & isbig(:, 3:end);
resptime = findRespTime(is3inrow);
figure;
hist(resptime, 0:0.1:2); %#ok<*HIST>
xlabel('Response time (sec)');
ylabel('Number of Units');
[sordid, srtidx] = sort(resptime);
figure;
imagesc(((-40:39)+0.5)*0.05, 1:numexcit, renormed(srtidx, :));
caxis([-3 3]);
colorbar;
axis xy;
xlabel('Time (sec)');
ylabel('Single units');
% Only include "good" response times, i.e. non-NaN, in the cumulative
% distribution function:
goodsordid = sordid(~isnan(sordid));
for k = 1:length(goodsordid)
    cumfunc(k) = k/length(goodsordid); %#ok<AGROW>
end
figure; plot(goodsordid, cumfunc);
xlabel('Response time (sec)');
ylabel(sprintf('Cumulative fraction\nof %d units', length(goodsordid)));

