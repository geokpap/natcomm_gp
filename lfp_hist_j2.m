function lfp_hist_j2(trials, avgflag, histflag, trialinfo)

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

    if avgflag || histflag ;
    end
    alltrialsamples = [];
    penis = 1;
    SaccTargetSums  = zeros(1,5);
    PeripheralTargetOrder = zeros(4,3);
    NumPeripheralTargsEntered = 0;
    GorBTargsEntered = 0;
    for trial = trials
        if lfp_SelectedTrials(trial) == 1 % do this loop only for selected trials
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
            
            % The time range of data we want to process is the union of the nominal
            % trial interval with the lfp_XLimAll interval, limited to the bounds
            % of the trial's recorded time segment.
            if isempty(lfp_XLimAll)
                xlimpoints = [0 0];
            else
                xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
            end
            if reftime ~= 0
                xlimabspoints = xlimpoints + refpoint;
            else
                % There was no reference event in this trial. The following becomes
                % a harmless redundant value when we compute startsample and
                % endsample:
                xlimabspoints = [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4)];
            end 
            startsample = max(...
                min( lfp_TrialIndex(trial,3), xlimabspoints(1)), ...
                lfp_TrialRec(trial,1) );
            endsample = min(...
                max( lfp_TrialIndex(trial,4), xlimabspoints(2)), ...
                lfp_TrialRec(trial,2) );
            trialsamplerange = startsample : endsample;
            if histflag 
                % Collect trial info to be used below for computing averaged
                % waveform. First find the sample index that is closest in time to
                % the lfp_AlignmentRef event:
                if (reftime - lfp_index2time(refpoint-1)) ...
                        < (lfp_index2time(refpoint) - reftime)
                    refpoint = refpoint - 1;
                end
                % Add to the trialinfo table, which contains the sample indices of
                % the start of trial, end of trial, and reference event:
                trialinfo = [ trialinfo
                    [ lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint ] ];
                
                % need to change these indeces for BDFormatNumber = 8
                trialstarttime = trialevents(find(trialevents(:,2) == 2),1);
                eyecalibration(5) = trialstarttime;
                
                % extract target position from BDOutput
                lfp_getTaskParams;
                
                % extract a range of sample indeces between two events of
                % interest, in order to plot the data in color (coding for time)
                % add option of user-input for these events
                coloroneventtime = trialevents(find(trialevents(:,2) == 16),1);
                coloroffeventtime = trialevents(find(trialevents(:,2) == 17),1);
                
                coloroneventindex = lfp_time2index(coloroneventtime);
                coloroffeventindex = lfp_time2index(coloroffeventtime);
                
                colortrialsamplerange = coloroneventindex : coloroffeventindex;
                greytrialsamplerange = startsample : (coloroneventindex - 1);
                blacktrialsamplerange = (coloroffeventindex + 1) : endsample;
                
                c = lfp_getEyeCalib(trial);
                calibrated_samples1c = ...
                    lfp_calibEye(lfp_Samples{1}(colortrialsamplerange), 'x', c);
                calibrated_samples2c = ...
                    lfp_calibEye(lfp_Samples{2}(colortrialsamplerange), 'y', c);
                
                % Minus sign because imagesc shows Y increasing down by
                % default:
                calibrated_samplesc = [calibrated_samples1c; -calibrated_samples2c]';
                
                % find maximum time range present in all trials
                pointsbefore = min(trialinfo(:,3) - trialinfo(:,1));
                pointsafter = min(trialinfo(:,2) - trialinfo(:,3));
                if ~isempty(lfp_XLimAll)
                    % convert lfp_XLimAll into points, truncate to maximum range
                    % present
                    pointsbefore = min(pointsbefore, ...
                        -round(lfp_XLimAll(1)/lfp_SamplePeriod));
                    pointsafter = min(pointsafter, ...
                        round(lfp_XLimAll(2)/lfp_SamplePeriod));
                end
                % compute averaged waveform for each filenum
                % avgwave = zeros(pointsbefore + pointsafter + 1, length(filenums));
                colnum = 0;
                
                allchannelsamples = [];
                
            end
            
            %distance from center of center target to center of peripheral
            %target
            ShapeSpacing = 200; % should read this value from BDoutput
            
            % extract target position from BDOutput
            lfp_getTaskParams2;
            
            Distance = [];
            for targetindex = 1:4
                
                switch (Target(targetindex))
                    case {0,1,3}
                        y_offset = 0;
                    case {5,6}
                        y_offset = -(2^0.5)*ShapeSpacing/2;
                    case {7,8}
                        y_offset = (2^0.5)*ShapeSpacing/2; %lower shapes should be negative
                    case {2}
                        y_offset = -ShapeSpacing;
                    case {4}
                        y_offset = 1*ShapeSpacing;
                end
                switch (Target(targetindex))
                    case {0,2,4}
                        x_offset = 0;
                    case {6,7}
                        x_offset = (2^0.5)*ShapeSpacing/2;
                    case {5,8}
                        x_offset = -(2^0.5)*ShapeSpacing/2;
                    case {3}
                        x_offset = ShapeSpacing;
                    case {1}
                        x_offset = -1*ShapeSpacing;
                end
                
                
                %calculate distance of samples from each of the four cues
                %(center, red, green, blue)
                
                Distance(targetindex,:) = (((calibrated_samplesc(:,1))' - x_offset).^2 + ((calibrated_samplesc(:,2))' - y_offset).^2).^0.5;
            end
             
            InTargets = Distance < 60;
            
            saccnum = 1;
            nextsample = 1;
            while ~isempty(InTargets)
                [targnum, samplenum] = find(InTargets); % DOUBLE-CHECK PT rad
                
                SaccTarget(saccnum) = targnum(1);
                
                nexttarg_index = find(targnum ~= SaccTarget(saccnum));
                if isempty(nexttarg_index)
                    InTargets = [];
                else
                    nextsample = samplenum(nexttarg_index(1));
                    InTargets = InTargets(:, nextsample : end);
                    saccnum = saccnum + 1;
                end
            end
            % SaccTarget
            
            % this matrix of zeroes and ones indicates whether each target
            % (Center, Red, Green, Blue) was entered at least once in this
            % trial
            boolean = zeros(1,3);
            for i = 2:4
                if ~isempty(find(SaccTarget == i))
                    boolean(i-1) = 1;
                end
            end
            
            NumPeripheralTargsEntered = sum(boolean) + NumPeripheralTargsEntered;
            
            % counts this trial if Green or Blue cues were entered
            if sum(boolean(2:3)) > 0
                GorBTargsEntered = 1 + GorBTargsEntered;
            end            
            penis = penis + 1;    
            
            % calculate #times each saccade target was entered in the
            % current trial            
            for targetindex = 1:4
                SaccTargetSums(targetindex) = sum(SaccTarget == targetindex) + SaccTargetSums(targetindex);
            end
            
            SaccTargetSums(targetindex+1) = length(SaccTarget) + SaccTargetSums(targetindex+1);
            
            PeripheralSaccTarget = SaccTarget(find(SaccTarget ~= 1));
            
            PeripheralTargetOrder(1,PeripheralSaccTarget(1)-1) = PeripheralTargetOrder(1,PeripheralSaccTarget(1)-1) + 1;
            PeripheralTargetOrder(4,PeripheralSaccTarget(end)-1) = PeripheralTargetOrder(4,PeripheralSaccTarget(end)-1) + 1;
            
            if length(PeripheralSaccTarget) > 2
                Interim_PeripheralSaccTarget = PeripheralSaccTarget(2:end-1); 
                i = 1;
                while  i <= length(Interim_PeripheralSaccTarget) && i <= 2
                    PeripheralTargetOrder(i+1,Interim_PeripheralSaccTarget(i)-1) = PeripheralTargetOrder(i+1,Interim_PeripheralSaccTarget(i)-1) + 1;
                    i = i + 1;
                end
            end         
        end    
    end
    AvgPeripheralTargetsEntered = NumPeripheralTargsEntered / penis
    PercentTrialsEnteredGorB = GorBTargsEntered / penis
    SaccTargetSums
    PeripheralTargetOrder
    
