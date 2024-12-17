
sessions=["2019-11-07_12-38-36"];
folders=dir('Z:\georgios\prez\ApAv');
%sessions={folders.name};
avgTables=containers.Map;
outcomeLabels={'AvoidApAv', 'ApproachApAv', 'ForcedRewardApAv', 'ForcedPunishApAv'};
for i=outcomeLabels
    avgTables(char(i))=table('Size',[length(sessions),4],'VariableTypes',["cell","cell","cell","cell"],'RowNames',sessions,'VariableNames',["ts","hr","hrv","medHRV"]);
end
eventTimeMatrix=[];
pavlovEventMatrix=[];
fullData=containers.Map;
itiData=containers.Map;
for s=sessions
    session=char(s);
    if ~contains(session,'20') %skip if its not a session
        continue
    elseif contains(session,'2018')
        continue
    end
    %check if this session was done in annex5
    %{
    if ~isempty(dir(['X:\georgios\Prez\Recordings\ApAv\' session 'fixed']))
            continue %if this session exists in annex5 don't redo it
    end
    if ~isempty(dir(['Z:\georgios\prez\ApAv\' session '\fixed\HR\'])) %check for a fixed folder
        path=['Z:\georgios\prez\ApAv\' session '\fixed\HR\'];
    else
        path=['Z:\georgios\prez\ApAv\' session '\HR\'];
    end
    %}
    path=['Z:\georgios\prez\ApAv\' session '\HR\'];
    %path=['X:\georgios\Prez\Recordings\ApAv\' session '\HR\'];
    [outcomes,averages,eventTimes,eventLabels,pavlovEvent,itis]=readSession(path);
    if isempty(outcomes)
        continue %if it was a bad session don't save its data
    end
    fullData(session)=outcomes;
    itiData(session)=itis;
    disp(eventTimes)
    disp(pavlovEvent)
    eventTimeMatrix=[eventTimeMatrix;eventTimes];
    pavlovEventMatrix=[pavlovEventMatrix;pavlovEvent];
    for i=outcomeLabels
        if ismember(i,keys(averages))
            tbl=avgTables(char(i));
            cat=averages(char(i));
            tbl(session,:)={{cat.TimeStamps},{cat.Heartrate},{cat.HRV},{cat.HRVMed}};
            avgTables(char(i))=tbl;
        end
    end
end

eventTimes=mean(eventTimeMatrix);
pavlovCueOff=mean(pavlovEventMatrix);

save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\prezAnnex4Plotting.mat',"avgTables","eventTimes","pavlovCueOff",'-v7.3')
save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\prezAnnex4Full.mat',"fullData","itiData",'-v7.3')


function [outcomes,averages,eventTimes,eventLabels,pavlovEvTime,itis]=readSession(pathnlx)
    msg=['reading ',pathnlx];
    disp(msg)
    %loading heartrate data
    %pathnlx='X:\georgios\LittleDebbie\Recordings\ApAvApAp\2022-02-07_10-02-40\HR\';
    matFileName='heart_rate_data.mat';
    pathToMatFile=fullfile(pathnlx,matFileName);
    load(pathToMatFile);
    
    %after reading in hr and rr check for nan values/empty
    if all(isnan(hr)) || all(isnan(rr))
        outcomes=[];
        averages=[];
        eventTimes=[];
        eventLabels=[];
        pavlovEvTime=[];
        itis=[];
        disp(['skipping ' pathnlx ' due to nan values'])
        return
    end

    %loading events data
    [ncs]=read_neuralynx_nev(pathnlx);
    
    %convert HRV from s to ms
    rrMS=rr*1000;
    
    %zero out timestamps and convert to seconds
    timeVec=[ncs.TimeStamp];
    zeroedTS=timeVec-timeVec(1);
    correctedEventsTS=double(zeroedTS)/1e+6;

    
    %seperate trials
    eventMarkers=[ncs.TTLValue];

    trialStarts=correctedEventsTS(eventMarkers==9); %find all events where event marker is 9
    
    %extract the times that correlate with key events
    [fixationOnsetTimes,fixationOnsetEvents]=getEventTimes(correctedEventsTS,eventMarkers,[214]);
    [cueOnsetTimes,cueOnsetEvents]=getEventTimes(correctedEventsTS,eventMarkers,[233 211 215]);
    [rewardOnsetTimes,rewardOnsetEvents]=getEventTimes(correctedEventsTS,eventMarkers,[244 239 238 217]);
    [inBetweenTimes,inBetweenEvents]=getEventTimes(correctedEventsTS,eventMarkers,[201,207,243,219,212,216,202,203,204,205,206,234]);
    [rwdLengthTimes,rwdLengthEvents]=getEventTimes(correctedEventsTS,eventMarkers,setdiff(linspace(1,200,200),[9 18]));
    
    roundedTS=round(ts,3); 
    ssevents=repmat(-1,1,length(roundedTS));
    for t=1:length(cueOnsetTimes)
        ssevents(find(roundedTS==round(cueOnsetTimes(t),3)))=cueOnsetEvents(t);%use find instead of logical bc we don't want to change size of ssevents
    end
    for t=1:length(rewardOnsetTimes)
        ssevents(find(roundedTS==round(rewardOnsetTimes(t),3)))=rewardOnsetEvents(t);
    end
    events=repmat(-1,1,length(roundedTS));
    for t=1:length(inBetweenTimes)
        events(find(roundedTS==round(inBetweenTimes(t),3)))=inBetweenEvents(t);
    end
    rewards=repmat(-1,1,length(roundedTS));
    for t=1:length(rwdLengthTimes)
        rewards(find(roundedTS==round(rwdLengthTimes(t),3)))=rwdLengthEvents(t);
    end
    endITIs=repmat(-1,1,length(roundedTS));
    for t=1:length(fixationOnsetEvents)
        endITIs(find(roundedTS==round(fixationOnsetTimes(t),3)))=fixationOnsetEvents(t);
    end
    %disp(size(ssevents))
    %create struct for each trial containing timestamp, event markers,heartrate and hrv
    %sampling rate of ts is approximatley every .001 second
    s=struct('TimeStamps',[],'Heartrate',[],'HRV',[],'FocusTimes',[],'Events',[],'RewardLengths',[],'EndITIs',[]);
    trials=repmat(s,1,length(trialStarts));
    
    for i=1:length(trialStarts)
        tStart=round(trialStarts(i),3); %round to 3dec bc ts only has 3 dec
        if i==length(trialStarts)
            tStop=round(ts(end),3);%there is no next trial so we use the end of ts to get the end of the trial
        else
            tStop=round(trialStarts(i+1),3); %start of next trial
        end
        %match start/stop with index in ts to get start
        idxStart= find(roundedTS==tStart);
        idxStop= find(roundedTS==tStop);
        if i~=length(trialStarts)
            idxStop=idxStop-1;
        end
        trials(i).TimeStamps=ts(idxStart:idxStop);%double check that duration isnt same for every trial
        trials(i).Heartrate=hr(idxStart:idxStop);%use that range to get hr and hrv
        trials(i).HRV=rrMS(idxStart:idxStop);
        trials(i).FocusTimes=ssevents(idxStart:idxStop);
        trials(i).Events=events(idxStart:idxStop);
        trials(i).RewardLengths=setdiff(rewards(idxStart:idxStop),-1,'stable');
        trials(i).EndITIs=endITIs(idxStart:idxStop);
    end
    
    %find timstamps where one of these markers occur, get this period for each trial
    %seperate based on trial type and further into outcome
    startingMarkers=containers.Map(['ApAv' 'RedForcedApAv' 'YellowForcedApAv'],[233 215 211]);
    outcomeMarkers=containers.Map({244, 239, 217, 238},{'AvoidApAv', 'ApproachApAv', 'ForcedRewardApAv', 'ForcedPunishApAv'});
    outcomes=containers.Map;
    itis=containers.Map;
    badTrials=0;
    for i=1:length(trials)
        trial=trials(i);
        %get events that aren't -1
        evs=find(trial.FocusTimes~=-1);
        if length(evs)<2
            badTrials=badTrials+1;
            continue
        elseif length(evs)>2
            timeDiffs=diff(evs);
            loc=find(timeDiffs==1);
            if loc==2
                evs=evs(1:2);
            elseif loc==1 | length(evs)==4
                evs=evs(2:3);
            else
                disp(evs)
                badTrials=badTrials+1;
                continue
            end 
        end
        %use their indeces to get the desired time frame
        startIdx=evs(1);
        stopIdx=evs(2);
        if ~ismember(trial.FocusTimes(startIdx),cell2mat(values(startingMarkers)))
            badTrials=badTrials+1;
            continue
        end
        if ~ismember(trial.FocusTimes(stopIdx),cell2mat(keys(outcomeMarkers)))
            badTrials=badTrials+1;
            continue
        end
        %get ITI data
        endITIIdx=find(trial.EndITIs~=-1);
        %index from begginging of data to the fixation onset
        ITIData=struct('TimeStamps',[],'Heartrate',[],'HRV',[],'TrialNum',[]);
        ITIData.TimeStamps=trial.TimeStamps(1:endITIIdx);
        ITIData.Heartrate=trial.Heartrate(1:endITIIdx);
        ITIData.HRV=trial.HRV(1:endITIIdx);
        ITIData.TrialNum=i;
        %get focus period data
        zoomedData=struct('TimeStamps',[],'Heartrate',[],'HRV',[],'Events',[],'Rewards',[],'TrialNum',[]);
        zoomedData.TimeStamps=trial.TimeStamps(startIdx:stopIdx);
        zoomedData.Heartrate=trial.Heartrate(startIdx:stopIdx);
        zoomedData.HRV=trial.HRV(startIdx:stopIdx);
        zoomedData.Events=trial.Events(startIdx:stopIdx);
        zoomedData.Rewards=trial.RewardLengths;
        zoomedData.TrialNum=i;
        %use then to classify and store the data
        category=outcomeMarkers(trial.FocusTimes(stopIdx));
        if ismember(category,keys(outcomes))
            old=outcomes(category);
            old(end+1)=zoomedData;
            outcomes(category)=old;

            oldITI=itis(category);
            oldITI(end+1)=ITIData;
            itis(category)=oldITI;
        else
            outcomes(category)=[zoomedData];
            itis(category)=[ITIData];
        end
    end
    
    %for each outcome, adjust all trials to be average trial length, then get average the data for each trial
    averages=containers.Map;
    categories=keys(outcomes);
    for o=1:numel(categories)
        cat=outcomes(char(categories(o)));
        %calculate average length of trials in the category for interpolation
        total=0;
        for i=1:length(cat)
            total=total+numel(cat(i).TimeStamps);
        end
        avTrialLen=round(total/length(cat));
        trialHrtrtMatrix=zeros(length(cat),avTrialLen);%create empty matrix where each row is a trial and each col is timepoint
        trialHRVMatrix=zeros(length(cat),avTrialLen);
        for i=1:length(cat)
            hrtrt=cat(i).Heartrate;
            cat(i).Heartrate=interp1(hrtrt,linspace(1,numel(hrtrt),avTrialLen));
            hrv=cat(i).HRV;
            cat(i).HRV=interp1(hrv,linspace(1,numel(hrv),avTrialLen));
            %add to category matrix
            trialHrtrtMatrix(i,:)=cat(i).Heartrate;
            trialHRVMatrix(i,:)=cat(i).HRV;
        end
        outcomes(char(categories(o)))=cat; %save interpolated data in outcomes
    
        %average over all trials for that category
        trialAvs=struct('TimeStamps',linspace(0,(avTrialLen-1)*.001,avTrialLen),'Heartrate',mean(trialHrtrtMatrix),'HRV',mean(trialHRVMatrix.'),'HRVMed',median(trialHRVMatrix));
        averages(char(categories(o)))=trialAvs;
    
    end
    
    %for each trial type need to find the average time the cross and square are displayed and the decision is made
    %ApAv
    try
        ApAvEvTimes=(findEvents(outcomes('ApproachApAv'),[234 219 201])+findEvents(outcomes('AvoidApAv'),[234 243 207]))/2;
        avTimeCSDisp=ApAvEvTimes(1); 
        avTimeDecision=ApAvEvTimes(2);
        avTimeHighlight=ApAvEvTimes(3);
        ApAvPavlovEvs=(findEvents(outcomes('ForcedRewardApAv'),[216])+findEvents(outcomes('ForcedPunishApAv'),[212]))/2;
    catch
        if outcomes.Count==0
            message=['no good trials in ',pathnlx];
            disp(message)
            outcomes=[];
            averages=[];
            eventTimes=[];
            eventLabels=[];
            pavlovEvTime=[];
            itis=[];
            return
        else
            message=['no ApAv trials in ',pathnlx];
            disp(message)
        end
    end

    eventTimes=[avTimeCSDisp,avTimeDecision,avTimeHighlight];
    pavlovEvTime=ApAvPavlovEvs;

    
    disp('bad trials: ')
    disp(badTrials)
    
    eventLabels={'Cross/Square Onset','Highlight On','Highlight Off'};

end

function evLocs=findEvents(trials,eventNums)
    evTotals=zeros(1,length(eventNums));
    evCounts=zeros(1,length(eventNums));
    for i=1:length(trials)
        trial=trials(i);
        for j=1:length(eventNums)
            ev=eventNums(j);
            idx=find(trial.Events==ev);
            tDiff=trial.TimeStamps(idx)-trial.TimeStamps(1);
            if length(idx)==1
                evTotals(j)=evTotals(j)+tDiff;
                evCounts(j)=evCounts(j)+1;
            end
        end
    end
    evLocs=(evTotals./evCounts);
end


function [times,events]=getEventTimes(ts,markers,targets)
    logical=ismember(markers,targets);
    times=ts(logical);
    events=markers(logical);
end
