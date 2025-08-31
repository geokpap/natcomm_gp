function test_dg_readEvents
% Regression test for dg_readEvents

%$Rev: 135 $
%$Date: 2011-12-01 19:44:37 -0500 (Thu, 01 Dec 2011) $
%$Author: dgibson $

load('dg_readEvents_testdata.mat');
evpath = which('s02acq02events.dat');
evpath2 = which('m47acq24events.nev');

[TStest, TTLtest, EStest, Hdrtest] = dg_readEvents(evpath);
if ~(isequal(TStest, TS) && isequal(TTLtest, TTL) ...
        && isequal(EStest, ES) && isequal(Hdrtest, Hdr))
    error('test_dg_readEvents:old', ...
        'Test failed on old file');
end
clear TStest TTLtest EStest Hdrtest

[TStest, TTLtest, EStest, Hdrtest] = dg_readEvents(evpath2);
if ~(isequal(TStest, TS2) && isequal(TTLtest, TTL2) ...
        && isequal(EStest, ES2) && isequal(Hdrtest, Hdr2))
    error('test_dg_readEvents:old', ...
        'Test failed on new file');
end

[TStest, TTLtest, EStest, Hdrtest] = dg_readEvents(evpath2, ...
    'mode', 2, [501 1000]);
if ~(isequal(TStest, TS2(501:1000)) && isequal(TTLtest, TTL2(501:1000)) ...
        && isequal(EStest, ES2(501:1000)) && isequal(Hdrtest, Hdr2))
    error('test_dg_readEvents:old', ...
        'Test failed for mode 2 on new file');
end

[TStest, TTLtest, EStest, Hdrtest] = dg_readEvents(evpath2, ...
    'mode', 4, [2130459200 3276504300]);
if ~(isequal(TStest, TS2(501:1000)) && isequal(TTLtest, TTL2(501:1000)) ...
        && isequal(EStest, ES2(501:1000)) && isequal(Hdrtest, Hdr2))
    error('test_dg_readEvents:old', ...
        'Test failed for mode 4 on new file');
end
disp('Tests passed.');
