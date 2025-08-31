function dg_realignVideo(sessiondir)
% Based on Katy's kt_realignVideo_script

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

lfp_declareGlobals;

%load session
evtsfound = false;
evtfilenames = {
    'events_90.evtsav' 
    'events.evtsav' 
    'events.nev' 
    'events.dat'
    };
for k = 1:length(evtfilenames)
    if exist(fullfile(sessiondir, evtfilenames{k}), 'file')
        lfp_read2('preset', sessiondir, {evtfilenames{k} ''});
        evtsfound = true;
        break
    end
end
if ~evtsfound
    error('dg_realignVideo:noevts', ...
        'Session %s contains no events file.', sessiondir);
    return
end

vtfilenames = {'vt1.nvt', 'vt1.dat'};
for k = 1:length(vtfilenames)
    infilename = fullfile(sessiondir, vtfilenames{k});
    if exist(infilename, 'file')
        break
    else
        infilename = '';
    end
end
if isempty(infilename)
    error('dg_realignVideo:nofile', ...
        '%s contains no video file', ...
        sessiondir );
end
[VT_timestamps, VT_X, VT_Y] = nlx2matvt_v3( ...
    infilename, [1 1 1 0 0 0], 0, 1, [] );
VT_timestamps = round(VT_timestamps / 100);

%for each trial, find time of Rec Off event, find closest break in
%VT timestamps, adjust VT timestamps
VTStartInd = 1;
NewVT_timestamps = [];
Offsets = [];
for trial = 1:length(lfp_SelectedTrials)
    if size(lfp_Events,1) < (lfp_TrialIndex(trial,2) + 3)
        TrialEvents = lfp_Events(lfp_TrialIndex(trial,1):end,2);
    else
        TrialEvents = lfp_Events(lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2)+3,2);
    end
    RecOffEventInd = lfp_TrialIndex(trial,1) + find(TrialEvents==2,1,'last') - 1;
    VTEndInd = find(VT_timestamps < (lfp_Events(RecOffEventInd,1)+5)*1e4, 1, 'last');
    TSoffset = round(lfp_Events(RecOffEventInd,1)*1e4) - VT_timestamps(VTEndInd);
    NewVT_timestamps(VTStartInd:VTEndInd) = VT_timestamps(VTStartInd:VTEndInd) + TSoffset;
    VTStartInd = VTEndInd+1;
    Offsets(trial) = TSoffset;
end

%convert to VACQ format
%open output file
[animaldir, sessionID ] = fileparts(sessiondir);
sessionID = upper(sessionID);
outfilename = fullfile(animaldir, ['V' sessionID '.DAT']);
if exist(outfilename, 'file')
    bakfilename = sprintf('%s.bak', outfilename);
    n = 0;
    while exist(bakfilename, 'file')
        n = n + 1;
        bakfilename = sprintf('%s.bak%d', outfilename, n);
    end
    movefile(outfilename, bakfilename);
end

fid = fopen(outfilename, 'w');
if fid == -1
    error('dg_realignVideo:output', ...
        'Could not open output file %s', ...
        outfilename );
end

try
    %write header
    fprintf(fid, '## %s\t%s\t (Created by MATLAB dg_realignVideo.m)\n', ...
        animaldir, sessionID);
    fprintf(fid, '##   Time Stamp\t\tX\t\tY\t\tFrame\t\tTrial\n');
    fprintf(fid, '##\n');
    %write video data
    for frame = 1:length(NewVT_timestamps)
        trialnum = 1;
        fprintf(fid,'\t%10d\t\t%d\t\t%d\t\t%d\t\t%d\n', ...
            round(NewVT_timestamps(frame)), ...
            VT_X(frame), VT_Y(frame), frame, trialnum );
    end
    fclose(fid);
catch
    fclose(fid);
end


