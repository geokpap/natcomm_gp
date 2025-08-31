function result = dg_parallellfplogstresstest(numcores, numreps)
for corenum = 1:numcores
    logname = 'junk0';
    k = 0;
    logfilename = [logname '.log'];
    while exist(logfilename, 'file')
        k = k + 1;
        logname = sprintf('junk%d', k);
        logfilename = [logname '.log'];
    end
    cmdstr = sprintf( ...
        'matlab -nosplash -nojvm -r "dg_lfplogstresstest(%d, sprintf(''%%06d'', dg_pid)); exit" -logfile %s & ', ...
        numreps, logfilename);
    result(corenum) = system(cmdstr);
end
