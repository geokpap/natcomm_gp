function test_dg_readCSC
% Regression test for dg_readCSC

%$Rev: 136 $
%$Date: 2011-12-06 17:06:17 -0500 (Tue, 06 Dec 2011) $
%$Author: dgibson $

filename{1} = 's02acq02lfp1.dat';
filename{2} = 'h24acq16csc1.ncs';

for fileidx = 1:2
    testpath = which(filename{fileidx});
    [p, n, e] = fileparts(testpath); %#ok<NASGU>
    matpath = fullfile(p, [n '.mat']);
    load(matpath);
    [TStest, Samplestest, Hdrtest] = dg_readCSC(testpath);
    % make insensitive to presence/absence of terminal empty cell in Hdr:
    if isempty(Hdrtest{end})
        Hdrtest(end) = [];
    end
    if ~(isequal(TStest, TS) && isequal(Samplestest, Samples) ...
            && isequal(Hdrtest, Hdr)) %#ok<NODEF>
        error('test_dg_readCSC:mode1', ...
            'Test failed on file %s', testpath);
    end
    clear TStest Samplestest Hdrtest
end

[TStest, Samplestest, Hdrtest] = dg_readCSC(testpath, ...
    'mode', 2, [501 1000]);
if ~(isequal(TStest, TS(501:1000)) && ...
        isequal(Samplestest, Samples(:,501:1000)))
    error('test_dg_readCSC:mode2', ...
        'Test failed for mode 2 on file %s', testpath);
end

start = [
    1143851587
    1143851588
    1143851589
    1143851600
    1143860000
    1143900000
    1144000000
    1144010000
    1144020000
    1144020500
    1144020548
    ];
for testnum = 1:length(start)
    timerange = [start(testnum) 1228500500];
    firstrecnum = find(TS >= timerange(1), 1);
    [TStest, Samplestest, Hdrtest] = dg_readCSC(testpath, ...
        'mode', 4, timerange);
    if ~(isequal(TStest, TS(firstrecnum:1000)) && ...
            isequal(Samplestest, Samples(:,firstrecnum:1000)))
        error('test_dg_readCSC:mode4', ...
            'Test failed for mode 4 on file %s', testpath);
    end
end

disp('Tests passed.');
