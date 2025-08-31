function [trials, numTrials, fileHeader, filename] = dj_readDataFile(varargin)
%   function [trials, numTrials, fileHeader, filename] = readDataFile(varargin)
%   varargin   -- optional variable.
%		  1st argument -- string, the file name of the data file.
%   trials     -- the information about all trials contained in the data file.
%   numTrials  -- total number of trials.
%   fileHeader -- file information.
%   filename   -- the file name from which the data is read in.
%   For details read the following note.
%
%   Note
%
%   This function reads the binary data supplied by Naotaka. The file
%   contains information about one trial. The file is binary with the following format:
%   (1) File header, consists of 8 words = 16 bytes. The meanings of the 8
%   words are: FFFF (marker), year, day, time, channal number, cluster,
%   corr, task id. These values will be read in and stored in the structure
%   fileHeader, which has the fields:
%       year,
%       day,
%       time,
%       taskID.
%   (2) Trial header. The length of the trial header is equal to 18 + the size
%   of the clusters. Here the size of the clusters is the sum of the 5th to 8th words in the file
%   header. The meaning of all the word are: EEEE (marker), 0, error (if 1,
%   error trial, 0 otherwise), number of spikes in cluster 1, number of
%   spikes in cluster 2, ..., number of spikes in the last cluster,
%   0,0,0,0,0,0, directions of saccades No. 1 & 2(lower 4 bits) + with or
%   without yellow (upper 4 bits), directions of saccades No. 3 & 4 (lower
%   4 bits) + target off timing (upper 4 bits), directions No. 5 & 6 (lower
%   4 bits) + interval information, task condition information, information
%   on number of saccades (lower 4 bits) + grid on/off (upper 4 bits).
%   (3) spike times of each cluster. Block for cluster 1, block for cluster
%   2, ..., block for the last cluster. The block for a cluster starts
%   withh DDD1(2,3,4,..), then time for each spike. note that clusters with
%   no spikes are also represented (the block has just the marker itself) .
%   (4) times of the events. starts with BBBB, 129 words.
%   (5) Repeat (2)-(4) until there is no more to read.
%   The information of a trial is strored in the array of structures
%   trials. Each trials have  the following fieds:
%       error (1 for error, 0 for correct, -100 for invalid numbers read in., -1 for event time not recorded properly),
%       numClusters,
%       numSpikes (in each clusters),
%       spikeTimes (1D cell array),
%       saccadeDirections (0-4. 0->up,1->right,2->down,3->left),
%       numSaccades,
%       yellow (1 yelllow on, 0 off),
%       targetOffCode,
%       targetOffTime,
%       interval,
%       condition,
%       gridOnOff,
%       taskStartTime,
%       fixationStartTime,
%       goStartTime (1D array, one element for each saccade),
%       saccadeOnsetTime (1D array),
%       saccadeOffsetTime (1D array),
%       rewardTime.
%
%   Written by Dezhe Jin, djin@mit.edu
%   Date 12/14/2002.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

% initialize.

numTrials = 0;
fileHeader = [];
filename = '';

% first get the data file.

if isempty(varargin)
    [filename, dirc] = uigetfile('/home/djin/dataDrive/Naotaka/*.dwd');
    if isempty(filename)
        disp('ERROR: No file is specified.');
        return;
    end
    filename = [dirc,filename];
else
    filename = varargin{1};
    disp(['Reading data from ',filename]);
end
fid = fopen(filename,'r');
if fid == -1
    trials = struct([]);
    disp(['ERROR: File could not be open. Check the file name.']);
    return;
end

% read the file header
fileHeaderArray = fread(fid,8,'uint16');
if  fileHeaderArray(1) ~= hex2dec('ffff')       % check the marker.
    trials = struct([]);
    disp('ERROR: Marker detection error. File header marker not seen.');
    return;
end
cluster = sum(fileHeaderArray(5:8));        	% total cluster number, only possible slots.
fileHeader = struct('year',fileHeaderArray(2),'day',fileHeaderArray(3),'time',fileHeaderArray(4),'taskID',fileHeaderArray(8));

% retrieve information about all trials.

numTrials = 0;
clear numSpikesClusterTrials;
while 1
    clear trialHeaderArray spikeArray saccadeDirections spikeTimes;
    validTrial = 1;

    % read the trial header
    trialHeaderArray = fread(fid,18+cluster,'uint16');

    if isempty(trialHeaderArray)
        break;
    end
    if trialHeaderArray(1) ~= hex2dec('eeee')
        disp('ERROR: Marker detection error. Trial header marker not seen.');
        error = -100;
        continue;
    end

    error = trialHeaderArray(3);

    number = trialHeaderArray(cluster+10);
    if number < 98
        validTrial = 0;
    else
        saccadeDirections(1:2) = getSaccadeCode(bitand(number-2,15));
        yellow = bitshift(number-98,-4);
    end

    number = trialHeaderArray(cluster+11);
    if number < 2
        validTrial = 0;
    else
        saccadeDirections(3:4) = getSaccadeCode(bitand(number-2,15));
        targetOffCode = bitshift(number,-4);
        if targetOffCode == 0 | targetOffCode == 1
            targetOffCode = targetOffCode + 2;
        elseif targetOffCode == 6 | targetOffCode == 7
            targetOffCode = targetOffCode -6;
        else
            targetOffCode = -1; % something is wrong.
        end
    end

    number = trialHeaderArray(cluster+12);
    if number < 2
        validTrial = 0;
    else
        saccadeDirections(5:6) = getSaccadeCode(bitand(number-2,15));
        interval0 = bitshift(number,-4);
        if  interval0== 6		% interval 0
            interval = 4444;
        elseif interval0 == 7	% interval 1
            interval = 4844;
        elseif interval0 == 0	% interval 2
            interval = 8444;
        elseif interval0 == 1	% interval 3
            interval = 6666;
        elseif interval0 == 2	% interval 4
            interval = 8888;
        end
    end
    condition = trialHeaderArray(cluster+13) - hex2dec('62');

    number = trialHeaderArray(cluster+14) - hex2dec('62') - trialHeaderArray(cluster+13) + hex2dec('80');
    if number < 0
        disp('WARNING: the number of saccades is not correctly read in. The trial is tagged error = -100.');
        validTrial = 0;
        numSaccades= -1;
        gridOnOff = -1;
    else
        numSaccades = bitand(number,15);
        gridOnOff = bitshift(number,-4);
    end

    % read the spike times of each block.
    spikeTimes = cell(1,cluster);
    for i=1:cluster
        spikeArray = fread(fid,trialHeaderArray(3+i)+1,'uint16');
        if spikeArray(1) == hex2dec('DDD0')+i
            spikeTimes(i) = {spikeArray(2:(trialHeaderArray(3+i)+1))'};
        else
            break;
        end
    end


    % read the events.
    eventArray = fread(fid,129,'uint16');
    if eventArray(1) ~= hex2dec('bbbb')
        disp('ERROR: Marker detection error. Event time marker not seen.');
        error = -100;
        continue;
    end
    taskStartTime = eventArray(2);
    fixationStartTime = eventArray(3);
    goStartTime = zeros(1,numSaccades);
    saccadeOnsetTime = zeros(1,numSaccades);
    saccadeOffsetTime = zeros(1,numSaccades);
    for i=1:numSaccades
        goStartTime(i) = eventArray(hex2dec('40')+i);
        saccadeOnsetTime(i) = eventArray(hex2dec('20')+2*(i-1)+1);
        saccadeOffsetTime(i) = eventArray(hex2dec('20')+2*(i-1)+2);
    end
    targetOffTime = eventArray(hex2dec('40')+numSaccades+2);
    if numSaccades > 0 & (targetOffTime > 10000 | targetOffTime <= 0)	% target off time is not recorded correctly. try recover.
        if targetOffCode == 0
            targetOffTime = saccadeOffsetTime(numSaccades);
        elseif targetOffCode == 1
            targetOffTime = goStartTime(numSaccades) + 400;
        end
    end

    rewardTime = eventArray(hex2dec('7f')+1);
    if (max(eventArray([2:3 (hex2dec('40')+1):(hex2dec('40')+numSaccades) ...
            (hex2dec('20')+1):(hex2dec('20')+2*numSaccades)])) > 10000 )
        error = -1;	% event times are not recorded correctly.
    end

    if validTrial == 1
        if numSaccades > 6
            saccadeDirections(7:8) = saccadeDirections(1:2);
        end
    else
        error = -100;		% mark invalid reading of the parameters.
        saccadeDirections = [];
        continue;
    end
    numTrials = numTrials + 1;
    trials(numTrials) = struct('error',error,'numClusters',cluster,...
        'numSpikes',trialHeaderArray(4:(3+cluster))',...
        'spikeTimes',{spikeTimes},'saccadeDirections',saccadeDirections,...
        'numSaccades',numSaccades,'yellow',yellow,...
        'targetOffCode',targetOffCode,'targetOffTime',targetOffTime,...
        'interval',interval,'condition',condition,'gridOnOff',gridOnOff,...
        'taskStartTime',taskStartTime,...
        'fixationStartTime',fixationStartTime,'goStartTime',goStartTime,...
        'saccadeOnsetTime',saccadeOnsetTime,...
        'saccadeOffsetTime',saccadeOffsetTime,'rewardTime',rewardTime);

    numSpikesClusterTrials(numTrials,:) = trialHeaderArray(4:13)';

end

% now determine the number of clusters and delete empty spikeTimes.

numValidClusters = max(find(sum(numSpikesClusterTrials) ~= 0));
for i=1:numTrials
    trials(i).numClusters = numValidClusters;
    trials(i).numSpikes   = numSpikesClusterTrials(i,1:numValidClusters);
    spikes = trials(i).spikeTimes(1:numValidClusters);
    trials(i).spikeTimes = spikes;
end

fclose(fid);




function [saccCode,saccCodeString]= getSaccadeCode(codeNum)
%	function saccCode = getSaccadeCode(codeNum)
%	This function converts the codeNum to saccade directions.
%	The function is mainly used in readDataFile.m
%	0 -> up, 1-> right, 2-> down, 3-> left.
%	Parameters:
%	codeNum  -- code number for the direction two consecutive saccades.
%	saccCode -- 1D array of two numbers for each directions.
%	saccCodeString -- letters for the directions.
%
%	Written by Dezhe Jin, djin@mit.edu
%	Date 1/6/2003.

if codeNum == 0
    saccCode = [0 1];
    saccCodeString = 'UR';
elseif codeNum == 1
    saccCode = [1 2];
    saccCodeString = 'RD';
elseif codeNum == 2
    saccCode = [2 3];
    saccCodeString = 'DL';
elseif codeNum == 3
    saccCode = [3 0];
    saccCodeString = 'LU';
elseif codeNum == 4
    saccCode = [1 0];
    saccCodeString = 'RU';
elseif codeNum == 5
    saccCode = [2 1];
    saccCodeString = 'DR';
elseif codeNum == 6
    saccCode = [3 2];
    saccCodeString = 'LD';
elseif codeNum == 7
    saccCode = [0 3];
    saccCodeString = 'UL';
else
    saccCode = [-1 -1];
    saccCodeString = 'WW';	% wrong code number!
end

