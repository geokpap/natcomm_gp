function lfp_saveVTfile(filename, VT_TS, VT_X, VT_Y)
% Saves Yasuo-format VT output file to <filename>.
% VT_TS - timestamps in seconds
% VT_X, VT_Y - coordinates in pixels

%$Rev:  $
%$Date:  $
%$Author:  $

lfp_declareGlobals;

%convert to VACQ format
fid = fopen(filename,'w+');
if fid == -1
    error('lfp_saveVTfile:badoutfile', ...
        'Could not open file %s for writing', ...
        filename );
end

try
    %write header
    fprintf(fid, '## %s (Created by MATLAB lfp_saveVTfile)\n', ...
        filename );
    fprintf(fid, '##   Time Stamp\t\tX\t\tY\t\tFrame\t\tTrial\n');
    fprintf(fid, '##\n');

    %write video data
    VT_TS = round(VT_TS * 1e4); % convert to Yasuo 0.1 ms units
    trialnum = 1;
    if size(lfp_TrialIndex,1) > 1
        % This is based on the assumption that the start of recording for a
        % trial is 2.0000 sec before its lfp_NominalTrialStart.  Note that
        % in real life this may result in a few orphaned frames at the
        % beginning or end of a trial that really belong to the adjacent
        % trial.
        nxtrialstart = round( ...
            lfp_Events(lfp_TrialIndex(2,1), 1) * 1e4 - 20000);
    else
        nxtrialstart = Inf;
    end
    for frame = 1:length(VT_TS)
        if trialnum < size(lfp_TrialIndex,1) && ...
                VT_TS(frame) >= nxtrialstart
            trialnum = trialnum + 1;
            if trialnum < size(lfp_TrialIndex,1)
                nxtrialstart = round( ...
                    lfp_Events(lfp_TrialIndex(trialnum,1), 1) * 1e4 - 20000);
            end
        end
        fprintf(fid,'\t%10d\t\t%d\t\t%d\t\t%d\t\t%d\n', ...
            VT_TS(frame), VT_X(frame), VT_Y(frame), ...
            frame, trialnum);
    end
catch
    s = lasterror;
    fclose(fid);
    rethrow(s);
end
fclose(fid);
