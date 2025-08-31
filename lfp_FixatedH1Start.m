function [H1StartFixations, a, b, c] = lfp_FixatedH1Start(trials, filenums, varargin)
% DAN: PLEASE ADD WINDOW OPTION!!
%H1StartFixations = lfp_FixatedH1Start(trials, filenums)
%H1StartFixations = lfp_FixatedH1Start(..., 'cal')
%H1StartFixations = lfp_FixatedH1Start(..., 'degsperpix', degsperpix)

% use: [H1StartFixations, numf5, numfc5, numcor5] = lfp_FixatedH1Start([], [1 2], 'cal');

% <trials>, <filenums> as in lfp_disp.
% Returns an
% array, H1StartFixations, containing information relating to the
% fixation surrounding H1Start event of each trial in lfp_SelectedTrials:
% whether there was a fixation as H1Start occurred,
% fixation start time, whether it was inside center cue, spatial location,
% distance from center, duration, "anticipation" duration, reaction time,
% CenterHold1 duration, trial# (not UniqueID!)
%
% Adapted from script "lfp_lib.77 pre-rel 2 JF\lfp_FixatedH1Start2.m"
%OPTIONS
%   'cal' - as in lfp_disp
%   'degsperpix' - sets scale factor to express positions in degrees;
%       default value is 1 so that position is in pixels.  Must be followed
%       by a number to give desired value.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

if strcmp(class(trials), 'char')
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

argnum = 1;
calflag = false;
degsperpix = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'cal'
            calflag = true;
        case 'degsperpix'
            degsperpix = varargin{argnum+1};
            argnum = argnum + 1;
        otherwise
            error('lfp_FixatedH1Start:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_FixatedH1Start:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
            num2str(filenums(find( ...
            ~ismember(filenums, lfp_ActiveFilenums) ))) ]);
end
if ~isequal(size(filenums), [1 2])
    error('lfp_FixatedH1Start:badFilenums', ...
        '<filenums> must be a 1x2 array' );
end
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_FixatedH1Start:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_FixatedH1Start:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials(find( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 ))) ]);
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_FixatedH1Start:badTrials2', '<trials> must be an integer vector.');
end

H1StartFixations = [];

FixatedH1Start = [];
FixatedCenter = [];
FixationDuration = [];
AnticipationDuration = [];

selectedtrial = 1; % just a counter for H1FixationStarts array

H1StartEvent = 16;
for trial = trials
    % do this loop for selected trials only
    if lfp_SelectedTrials(trial) == 1
        
        % extract target position from BDOutput
        lfp_getTaskParams;
        
        eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
        
        % find 'reftime', the lfp_AlignmentRef absolute timestamp:
        trialevents = lfp_Events(eventrange,:);
        reftime = trialevents( ...
            find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
            1 );
        if length(reftime) == 0
            reftime = 0;
        else
            reftime = reftime(1);
        end
        refpoint = lfp_time2index(reftime);
        
        if trialevents((find(trialevents(:,2) == H1StartEvent)) - 1, 2) == 132
            H1StartFixations(selectedtrial, 1) = 1; % eyes were fixated at H1Start
            
            % Collect trial info to be used below for computing averaged
            % waveform. First find the sample index that is closest in time to
            % the lfp_AlignmentRef event:
            if (reftime - lfp_index2time(refpoint-1)) ...
                    < (lfp_index2time(refpoint) - reftime)
                refpoint = refpoint - 1;
            end
            
            FixatedH1StartTime = trialevents(find(trialevents(:,2) == ...
                H1StartEvent),1);
            
            FixatedH1StartIndex = lfp_time2index(FixatedH1StartTime);
            
            if calflag
                calibrated_samples1c = degsperpix * ...
                    lfp_Samples{filenums(1)}(FixatedH1StartIndex);
                calibrated_samples2c = degsperpix * ...
                    lfp_Samples{filenums(2)}(FixatedH1StartIndex);
            else
                c = lfp_getEyeCalib(trial);
                calibrated_samples1c = degsperpix * ...
                    lfp_calibEye(lfp_Samples{filenums(1)}(FixatedH1StartIndex), 'x', c);
                calibrated_samples2c = degsperpix * ...
                    lfp_calibEye(lfp_Samples{filenums(2)}(FixatedH1StartIndex), 'y', c);
            end
            
            % short (0) or long (1) CenterHold1 delay
            H1StartFixations(selectedtrial, 9) = CenterHold1;
            
            H1StartFixations(selectedtrial, 3) = calibrated_samples1c;
            H1StartFixations(selectedtrial, 4) = calibrated_samples2c;
            
            Distance = (((calibrated_samples1c).^2 + ...
                (calibrated_samples2c).^2).^0.5);
            
            H1StartFixations(selectedtrial, 5) = Distance;
            
            % AnticipationDuration
            H1StartFixations(selectedtrial, 7) = FixatedH1StartTime - ...
                trialevents((find(trialevents(:,2) == H1StartEvent)) - 1, 1);
            
            % FixationDuration
            % change this to find the first fixationoff event post-H1Start
            H1StartFixations(selectedtrial, 6) = ...
                trialevents((find(trialevents(:,2) == H1StartEvent)) + 1, 1) - ...
                trialevents((find(trialevents(:,2) == H1StartEvent)) - 1, 1);
            
            % Are the eyes fixated within center cue?
            if ( Distance < 77)
                H1StartFixations(selectedtrial, 2) = 1; % eyes were fixated at CENTER at H1Start
            end
        H1StartFixations(selectedtrial, 10) = trial;    % Dan please make this UNIQUE ID
        selectedtrial = selectedtrial + 1;    
        end
        
    end
end

% RT = (FixationDuration - AnticipationDuration)
H1StartFixations(:, 8) = H1StartFixations(:, 6) - H1StartFixations(:, 7);

% mean RT of all FixatedH1Start trials
mean(H1StartFixations(find(H1StartFixations(:, 8)),8))

% find all NON-CENTER FixatedH1Start trials
% find(H1StartFixations(:,1) - H1StartFixations(:,2))

% mean RT of all NON-CENTER FixatedH1Start trials
mean(H1StartFixations(find(H1StartFixations(:,1) - H1StartFixations(:,2)), 8))

% mean RT of all CENTER FixatedH1Start trials
mean(H1StartFixations(find(H1StartFixations(:, 2)),8))

% % mean RT of all CENTER FixatedH1Start SHORT H1 trials
% mean(H1StartFixations(find((H1StartFixations(:,2) - H1StartFixations(:,9)) > 0), 8))

% mean RT of all CENTER FixatedH1Start SHORT H1 trials????


a = length(find(H1StartFixations(:,1)))

b = length(find(H1StartFixations(:,2)))

c = sum(lfp_SelectedTrials)