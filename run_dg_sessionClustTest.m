function run_dg_sessionClustTest(sessiondir)
% Runs 'dg_sessionClustTest' with the usual value of <retrig>, displays the
% text output, and then a talley of the number of clusters in each
% category.



retrig = 219e-6; % usual Georgios value
exceptions = {}; % usual Georgios value
[sesreport, sesreportary] = dg_sessionClustTest(sessiondir, ...
    retrig, exceptions);
disp(sesreport);
catvals = unique(sesreportary(:,15));
for k=1:length(catvals)
    slct=ismember(sesreportary(:,15), catvals{k});
    fprintf('%s: %d\n', catvals{k}, sum(slct));
end
