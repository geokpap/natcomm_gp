function test_lfp_gatherTrialData(addresult)
% Regression test for lfp_gatherTrialData.
%INPUTS 
% addresult: if true, compare to lfp_disp results and save the good results
%   to test_lfp_gathertrialdata.mat.


% To add more test cases, use getdispdata or getdispdatamultichannel to
% regression test against lfp_disp.  Note that the results can be off by
% one sample on either end and still be considered correct.
%   Helpful one-liners:
% isequal(myydata, data)
% isequal(myydata, data(:,2:end))
% refresult(dataidx, xlimidx, fnidx, tridx, winidx, 1:2) = struct('data', data, 'interval', interval, 'win', win)
% save('/Volumes/smbshare/software/sandbox/dgibson/dg_lib/test_lfp_gathertrialdata.mat', 'refresult');
% [dataidx, xlimidx, fnidx, tridx, winidx, optsidx]
% winidx=winidx+1; [ data, interval, win ] = lfp_gatherTrialData( trvals{tridx}, fnvals{fnidx}, winvals{winidx}, optsvals{optsidx}{:} );
% getdispdatamultichannel; isequal(myydata, data)
% isequal(myydata, data(:,2:end,:))

% See dg_lib for version history prior to 29-Aug-2013.


%$Rev: 305 $
%$Date: 2013-09-05 18:24:19 -0400 (Thu, 05 Sep 2013) $
%$Author: dgibson $

global lfp_AlignmentRef lfp_XLimAll lfp_Samples

% dataset 1 has a good assortment of different length trials, but it is all
% in on recording segment.  dataset 2 is nice and old (April 2002) but not
% SO old as to have no outcome codes (a problem with s02acq02), and has 93
% (!) recording segments.
datasets = { {'h24acq16.nev', 'h24acq16csc1.ncs', 'h24acq16csc2.ncs'}
    {'s17acq10events.dat' 's17acq10lfp1.dat' 's17acq10lfp2.dat'} };
xlimvals = {[] [-20 -10] [-2 -1] [-1 1] [1 2] [10 20]};
fnvals = {1 1:3 (1:3)'};
trvals = {10 11:15 (11:15)'};
winvals = xlimvals;
%optsvals = { {} {'pad'} {'logical'} {'pad', 'logical'} };
optsvals = { {} {'logical'} };

mypath = which('test_lfp_gatherTrialData');
mydir = fileparts(mypath);
lfp_changeSetup('katy');
load(fullfile(mydir, 'test_lfp_gathertrialdata.mat'));

for dataidx = 1:length(datasets)
    lfp_read2('preset', mydir, datasets{dataidx});
    lfp_createWave(@lfp_waverecord, 1, ...
        lfp_Samples{1}(:)>median(lfp_Samples{1}(:)));
    lfp_Samples{length(lfp_Samples)} = ...
        logical(lfp_Samples{length(lfp_Samples)});
    lfp_AlignmentRef = 14;
    for xlimidx = 1:length(xlimvals)
        for fnidx = 1:length(fnvals)
            for tridx = 1:length(trvals)
                for winidx = 1:length(winvals)
                    for optsidx = 1:length(optsvals)
                        lfp_XLimAll = xlimvals{xlimidx};
                        [ data, interval, win ] = ...
                            lfp_gatherTrialData( trvals{tridx}, ...
                            fnvals{fnidx}, winvals{winidx}, ...
                            optsvals{optsidx}{:} );
                        if all([dataidx, xlimidx, fnidx, tridx, ...
                                winidx, optsidx] <= size(refresult))
                            r = refresult(dataidx, xlimidx, fnidx, tridx, ...
                                winidx, optsidx);
                            if isequalwithequalnans(data, r.data) && ...
                                    isequalwithequalnans(interval, r.interval) ...
                                    && isequalwithequalnans(win, r.win)
                                fprintf('passed test: %s\n', dg_thing2str( ...
                                    [dataidx, xlimidx, fnidx, ...
                                    tridx, winidx, optsidx] ));
                            else
                                if addresult
                                    regressiontest_v_lfp_disp;
                                else
                                    error('test_lfp_gatherTrialData:bad', ...
                                        'Failed test [%d,%d,%d,%d,%d,%d].', ...
                                        dataidx, xlimidx, fnidx, tridx, winidx, ...
                                        optsidx);
                                end
                            end
                        else
                            % ran out of reference results
                            if addresult
                                regressiontest_v_lfp_disp;
                            else
                                error('Had to stop here.');
                            end
                        end
                    end
                end
            end
        end
    end
end
disp('Tests passed.');
if addresult
    save(fullfile(mydir, 'test_lfp_gathertrialdata.mat'), 'refresult');
end

    function regressiontest_v_lfp_disp
        myydata = [];
        if length(fnvals{fnidx}) > 1
            getdispdatamultichannel;
        else
            getdispdata;
        end
        if isempty(myydata) && isempty(data) || isequal(myydata, data) || ...
                ~isempty(data) && isequal(myydata, data(:,2:end,:)) || ...
                ~isempty(myydata) && isequal(data(1,:,:), myydata(1,2:end,:))
            refresult(dataidx, xlimidx, fnidx, ...
                tridx, winidx, optsidx) = ...
                struct('data', data, ...
                'interval', interval, 'win', win);
            fprintf('added refresult(%s)\n', dg_thing2str( ...
                [dataidx, xlimidx, fnidx, ...
                tridx, winidx, optsidx] ));
        else
            save(fullfile(mydir, 'test_lfp_gathertrialdata.mat'), ...
                'refresult');
            error('no, this one really fails!');
        end
        
        
        function getdispdata
            % for single channel
            try
                hF = lfp_disp(trvals{tridx}, fnvals{fnidx}, ...
                    winvals{winidx}, 'ovr');
            catch e
                switch e.identifier
                    case 'lfp_disp:nodata2'
                        myydata = [];
                    otherwise
                        myydata = NaN;
                end
                return
            end
            hL = findobj(hF, 'Type', 'line');
            traceidx = cellfun(@length, get(hL, 'XData')) > 2;
            hL(~traceidx) = [];
            if ismember( 'logical', optsvals{optsidx})
                myydata = logical(get(hL(1), 'YData'));
            else
                myydata = get(hL(1), 'YData');
            end
            for k = 2:length(hL)
                if ismember( 'logical', optsvals{optsidx})
                    myydata(k,:) = get(hL(k), 'YData');
                else
                    myydata(k,:) = get(hL(k), 'YData');
                end
            end
            myydata = flipud(myydata);
            close(hF);
        end
        
        function getdispdatamultichannel
            try
                hF = lfp_disp(trvals{tridx}, fnvals{fnidx}, ...
                    winvals{winidx}, 'ovr');
            catch e
                switch e.identifier
                    case 'lfp_disp:nodata2'
                        myydata = [];
                    otherwise
                        myydata = NaN;
                end
                return
            end
            hA = findobj(hF, 'Type', 'axes');
            for m = 1:length(hA)
                hL = findobj(hA(m), 'Type', 'line');
                traceidx = cellfun(@length, get(hL, 'XData')) > 2;
                hL(~traceidx) = [];
                if isempty(myydata)
                    if ismember( 'logical', optsvals{optsidx})
                        myydata = logical(get(hL(1), 'YData'));
                    else
                        myydata = get(hL(1), 'YData');
                    end
                end
                % filenums X timepoints X trials format
                for k = 1:length(hL)
                    if ismember( 'logical', optsvals{optsidx})
                        myydata(m,:,k) = logical(get(hL(k), 'YData'));
                    else
                        myydata(m,:,k) = get(hL(k), 'YData');
                    end
                end
            end
            myydata = flipdim(myydata,1);
            myydata = flipdim(myydata,3);
            close(hF);
        end
    end

end
