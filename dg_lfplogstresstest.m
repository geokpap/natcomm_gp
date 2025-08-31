function dg_lfplogstresstest(nmsgs, idstr)
% To thoroughly exercise lfp_log.

%$Rev: 81 $
%$Date: 2010-09-20 17:13:28 -0400 (Mon, 20 Sep 2010) $
%$Author: dgibson $

global lfp_LogFileName
if isempty(lfp_LogFileName)
    logname = 'lfp_lib.log';
    lfp_LogFileName = which(logname);
    if isempty(lfp_LogFileName)
        lfp_LogFileName = fullfile(pwd, logname);
    end
end
for k = 1:nmsgs
    lfp_log(sprintf( ...
        'idstr: %s #%d for all good men to come to the aid of their gerbil.', ...
        idstr, k));
    msg = sprintf('idstr: %s #%d some random lines of text:', idstr, k');
    for j = 1:10
        msg = sprintf('%s\n%s', msg, char(round(94*rand(1,60)+32)));
    end
    lfp_log(msg);
end
