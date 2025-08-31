function [lfp_hhtOneTrial_result, lfp_hhtOneTrial_fedge, ...
        lfp_hhtOneTrial_f, lfp_hhtOneTrial_A, lfp_hhtOneTrial_imf] = ...
        lfp_hhtOneTrial(idstr)
% Ridiculous kludge to parallelize at the level of trials without moving to
% Matlab Distributed Computing Server.
%   If <idstr> is a string, then this function runs in batch mode, reading
% input from and writing output to files.
%   If <idstr> is a cell array, then it should be a vector containing
% values for the variables 'waveform', 'emd_opts', 'hht_opts', 'dt',
% 'plotfreqsflag', 'freqedges' in that order, and the output is returned in
% the normal way.  If idstr{1} (the value of 'waveform') is a cell array,
% then its first element is assumed to contain one set of IMFs as returned
% in each cell of <imfs> from lfp_hht, and they are used instead of calling
% dg_emd.

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

global lfp_LogFileName

if isempty(lfp_LogFileName)
    logname = 'lfp_lib.log';
    lfp_LogFileName = which(logname);
    if isempty(lfp_LogFileName)
        lfp_LogFileName = fullfile(pwd, logname);
    end
end
if ischar(idstr)
    msg = sprintf('lfp_hhtOneTrial starting pid %d; idstr:%s', dg_pid, idstr);
    disp(msg);
    lfp_log(msg);
    inputfilename = sprintf('%s_input.mat', idstr);
    load(inputfilename, '-mat');
else
    waveform = idstr{1};
    emd_opts = idstr{2};
    hht_opts = idstr{3};
    dt = idstr{4};
    plotfreqsflag = idstr{5};
    freqedges = idstr{6};
end

lfp_hhtOneTrial_fedge = [];
lfp_hhtOneTrial_f = [];
lfp_hhtOneTrial_A = [];
if iscell(waveform)
    lfp_hhtOneTrial_imf = waveform{1};
else
    lfp_hhtOneTrial_imf = dg_emd(waveform, emd_opts{:});
end
lfp_hhtOneTrial_imf(end,:) = [];    % the last 'component' is residual junk
try
    if ischar(idstr)
        fprintf(1, '%s %s\n', datestr(now, 0), 'Running dg_hht');
    end
    dg_hhtbadfreqs = warning('off', 'dg_hht:badfreqs');
    [S, lfp_hhtOneTrial_fedge, lfp_hhtOneTrial_f, ...
        lfp_hhtOneTrial_A] = dg_hht(lfp_hhtOneTrial_imf, ...
        'freqedges', freqedges, hht_opts{:});
    warning(dg_hhtbadfreqs.state, 'dg_hht:badfreqs');
    if ischar(idstr)
        fprintf(1, '%s %s\n', datestr(now, 0), 'Done dg_hht');
    end
    % Remove the initial NaNs:
    S(:,1) = [];
    % Condense each group of <dt> columns to one column by summing:
    if dt > 1
        % Each column in the aggregate will contain columns (n*dt) + (1:dt)
        % for n ranging from 0 to {number of aggregate columns} - 1.  So
        % the very last column going into the aggregate will be an integral
        % multiple of dt.
        lastpt = floor(size(S,2) / dt) * dt;
    else
        lastpt = size(S,2);
    end
    lfp_hhtOneTrial_result = S(:,1:dt:lastpt);
    for k = 2:dt
        lfp_hhtOneTrial_result = lfp_hhtOneTrial_result + S(:,k:dt:lastpt);
    end
    lfp_hhtOneTrial_result = lfp_hhtOneTrial_result/dt;
    if ischar(idstr)
        fprintf(1, '%s %s\n', datestr(now, 0), 'Saving result');
        outputfilename = sprintf('%s_output.mat', idstr);
        save(outputfilename, 'lfp_hhtOneTrial_result', ...
            'lfp_hhtOneTrial_fedge', ...
            'lfp_hhtOneTrial_f', 'lfp_hhtOneTrial_A', ...
            'lfp_hhtOneTrial_imf');
    end
catch s
    if ischar(idstr)
        fprintf(1, '%s %s\n', datestr(now, 0), 'Caught error');
        save([idstr '_err.mat'], 's');
        exit;
    else
        rethrow(s);
    end
end
if ischar(idstr)
    msg = sprintf('lfp_hhtOneTrial pid %d Setting done flag %s', ...
        dg_pid, datestr(now, 0));
    disp(msg);
    lfp_log(msg);
    fid = fopen(sprintf('%s.done', idstr), 'w');
    fclose(fid);
end

