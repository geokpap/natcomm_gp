function [counts, avgdist, density] = dg_showClustStats(idx, C, sumD)
%   <idx>, <C> are as returned by dg_seqKMeans.
%   <counts>:  number of points in each cluster
%   <avgdist>:  total distance divided by number of points for each
%       cluster
%   <density>:  number of points divided by avgdist^3 (Note that the term
%       "density" is used loosely here, as it only applies literally to
%       spherical clusters)

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

for k=1:length(C)
    counts(k) = sum(idx==k);
    avgdist(k) = sumD(k)/sum(idx==k);
end
density = counts ./ (avgdist .^ 3);

figure;
subplot(2,2,1); 
plot(counts, ':*');
xlabel('cluster #');
ylabel('counts');
subplot(2,2,2);
plot(avgdist, ':*');
xlabel('cluster #');
ylabel('dist/memb');
subplot(2,2,3);
plot(density, ':*');
xlabel('cluster #');
ylabel('density');
subplot(2,2,4);
hist(density);
xlabel('density');
ylabel('clusters');
