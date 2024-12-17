%sessions=["2020-03-21_concatfixed","2021-12-17_concat","2021-12-20_concat","2022-01-06_concat","2022-05-10_concatfixed","2022-06-07_concat","2022-07-21_concat","2022-08-26_concat"];
%sessions=["2021-12-17_concat"];
folders=dir('X:\georgios\LittleDebbie\Recordings\ApAv');
sessions={folders.name};
avgTables=containers.Map;
outcomeLabels={'AvoidApAv', 'ApproachApAv', 'ForcedRewardApAv', 'ForcedPunishApAv'};
for i=outcomeLabels
    avgTables(char(i))=table('Size',[length(sessions),2],'VariableTypes',["cell","cell"],'RowNames',sessions,'VariableNames',["ts","lick"]);
end
eventTimeMatrix=[];
pavlovEventMatrix=[];
zoomedData=containers.Map;
fullData=containers.Map;
for s=sessions
    %continue
    session=char(s);
    if ~contains(session,'20') %skip if its not a session
        continue
    elseif contains(session,'2018')
        continue
    %elseif contains(session,'2019-11-07_12-38-36')
        %continue
    end
    %{
    %check if this session was done in annex5
    if ~isempty(dir(['X:\georgios\Prez\Recordings\ApAv\' session 'fixed']))
            continue %if this session exists in annex5 don't redo it
    end
    if ~isempty(dir(['Z:\georgios\prez\ApAv\' session '\fixed\Lick\'])) %check for a fixed folder
        path=['Z:\georgios\prez\ApAv\' session '\fixed\Lick\'];
    else
        path=['Z:\georgios\prez\ApAv\' session '\Lick\'];
    end
    %}
    path=['X:\georgios\LittleDebbie\Recordings\ApAv\' session '\lick\'];
    [outcomes,averages,eventTimes,eventLabels,pavlovEvent,trials]=readSession(path);
    if isempty(outcomes)
        continue %if it was a bad session don't save its data
    end
    zoomedData(session)=outcomes;
    fullData(session)=trials;
    disp(eventTimes)
    disp(pavlovEvent)
    eventTimeMatrix=[eventTimeMatrix;eventTimes];
    pavlovEventMatrix=[pavlovEventMatrix;pavlovEvent];
    for i=outcomeLabels
        if ismember(i,keys(averages))
            tbl=avgTables(char(i));
            cat=averages(char(i));
            tbl(session,:)={{cat.TimeStamps},{cat.Lick}};
            avgTables(char(i))=tbl;
        end
    end
end

%readSession('X:\georgios\Prez\Recordings\ApAv\2019-07-16_15-08-48fixed\Lick\')

eventTimes=mean(eventTimeMatrix);
pavlovCueOff=mean(pavlovEventMatrix);

%save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\debbieConcatPlotting.mat',"avgTables","eventTimes","pavlovCueOff",'-v7.3')
%save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\debbieConcatZoomed.mat',"zoomedData",'-v7.3')
%save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\debbieConcatFull.mat',"fullData",'-v7.3')

save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\prezAnnex4Plotting.mat',"avgTables","eventTimes","pavlovCueOff",'-v7.3')
save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\prezAnnex4Zoomed.mat',"zoomedData",'-v7.3')
save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\prezAnnex4Full.mat',"fullData",'-v7.3')

function [outcomes,averages,eventTimes,eventLabels,pavlovEvTime,trials]=readSession(pathnlx)
    msg=['reading ',pathnlx];
    disp(msg)
    %load licking data
    %pathnlx='X:\georgios\LittleDebbie\Recordings\ApAv\2022-08-09_07-58-56\lick';
    matFileName='lick_data.mat';
    pathToMatFile=fullfile(pathnlx,matFileName);
    try
        load(pathToMatFile);
    catch
        outcomes=[];
        averages=[];
        eventTimes=[];
        eventLabels=[];
        pavlovEvTime=[];
        trials=[];
        disp(['skipping ' pathnlx ' due to no lick file'])
        return
    end
    
    %after reading in hr and rr check for nan values/empty
    if all(isnan(samples))
        outcomes=[];
        averages=[];
        eventTimes=[];
        eventLabels=[];
        pavlovEvTime=[];
        trials=[];
        disp(['skipping ' pathnlx ' due to nan values'])
        return
    end

    %loading events data
    [ncs]=read_neuralynx_nev(pathnlx);
    
    %zero out timestamps and convert to seconds
    timeVec=[ncs.TimeStamp];
    zeroedTS=timeVec-timeVec(1);
    correctedEventsTS=double(zeroedTS)/1e+6;
    
    DIFF=ts(2)-ts(1);

    eventMarkers=[ncs.TTLValue];
    cueOnsetEvents=[233 211 215];
    rewardOnsetEvents=[244 239 238 217];
    ssevents=repmat(-1,1,length(ts));
    events=repmat(-1,1,length(ts));
    if contains(pathnlx,'concat') %find index where time was reset so it can be corrected
        diffs=diff(ts);
        RESTART=find(diffs>1);
    end
    
    first=false;
    for i=1:length(eventMarkers)
        time=correctedEventsTS(i);
        if contains(pathnlx,'concat') && time>RESTART*DIFF
            if first==false
                idx=RESTART+1;
                first=true;
            else
                idx=int64((time-(ts(RESTART+1)-ts(RESTART)))/DIFF);
            end
        else
            idx=int64(time/DIFF);
        end
        ev=eventMarkers(i);    
        if ismember(ev,cueOnsetEvents) | ismember(ev,rewardOnsetEvents)
            ssevents(idx+1)=ev;
        else
            events(idx+1)=ev;
        end
        if length(events)~=length(ts) %4297601-4297600
            disp('what')
        end
    end
    
    
    %create struct for each trial containing timestamp, event markers, and lick data
    %sampling rate of ts is approximatley every .001 second
    correctedTS=ts-ts(1);
    trialStarts=find(events==9); %find all events where event marker is 9
    s=struct('TimeStamps',[],'Lick',[],'FocusTimes',[],'Events',[],'RewardLengths',[],"TrialNum",[]);
    trials=repmat(s,1,length(trialStarts));
    for i=1:length(trialStarts)
        idxStart=trialStarts(i); %round to 3dec bc ts only has 3 dec
        if i==length(trialStarts)
            idxStop=length(correctedTS);%there is no next trial so we use the end of ts to get the end of the trial
        else
            idxStop=trialStarts(i+1); %start of next trial
        end
        if i~=length(trialStarts)
            idxStop=idxStop-1;
        end
        %get rewards
        trialTime=correctedTS(idxStart:idxStop);
        trialEvents=events(idxStart:idxStop);
        rewards=logical((trialEvents<200).*(trialEvents~=-1).*(trialEvents~=9).*(trialEvents~=18));
        trials(i).TimeStamps=trialTime;%double check that duration isnt same for every trial
        trials(i).Lick=samples(idxStart:idxStop);%use that range to get licking data
        trials(i).FocusTimes=ssevents(idxStart:idxStop);
        trials(i).Events=events(idxStart:idxStop);
        trials(i).RewardLengths=trialEvents(rewards);
        trials(i).TrialNum=i;
    end
    
    %find timstamps where one of these markers occur, get this period for each trial
    %seperate based on trial type and further into outcome
    startingMarkers=containers.Map(['ApAv' 'RedForcedApAv' 'YellowForcedApAv'],[233 215 211]);
    outcomeMarkers=containers.Map({244, 239, 217, 238},{'AvoidApAv', 'ApproachApAv', 'ForcedRewardApAv', 'ForcedPunishApAv'});
    outcomes=containers.Map;
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
                disp(trial.FocusTimes(evs))
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
        %get focus period data
        zoomedData=struct('TimeStamps',[],'Lick',[],'Events',[],'Rewards',[],'TrialNum',[]);
        zoomedData.TimeStamps=trial.TimeStamps(startIdx:stopIdx);
        zoomedData.Lick=trial.Lick(startIdx:stopIdx);
        zoomedData.Events=trial.Events(startIdx:stopIdx);
        zoomedData.Rewards=trial.RewardLengths;
        zoomedData.TrialNum=trial.TrialNum;
        if zoomedData.TrialNum==335
            disp('hi')
        end
        %use then to classify and store the data
        category=outcomeMarkers(trial.FocusTimes(stopIdx));
        if ismember(category,keys(outcomes))
            old=outcomes(category);
            old(end+1)=zoomedData;
            outcomes(category)=old;
        else
            outcomes(category)=[zoomedData];
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
        trialLickMatrix=zeros(length(cat),avTrialLen);%create empty matrix where each row is a trial and each col is timepoint
        for i=1:length(cat)
            l=cat(i).Lick;
            cat(i).Lick=interp1(l,linspace(1,numel(l),avTrialLen));
            %add to category matrix
            trialLickMatrix(i,:)=cat(i).Lick;
        end
        outcomes(char(categories(o)))=cat; %save interpolated data in outcomes
    
        %average over all trials for that category
        trialAvs=struct('TimeStamps',linspace(0,(avTrialLen-1)*.001,avTrialLen),'Lick',mean(trialLickMatrix));
        averages(char(categories(o)))=trialAvs;
    end
    
    ApAvEvTimes=(findEvents(outcomes('ApproachApAv'),[234 219 201])+findEvents(outcomes('AvoidApAv'),[234 243 207]))/2;
    avTimeCSDisp=ApAvEvTimes(1); 
    avTimeDecision=ApAvEvTimes(2);
    avTimeHighlight=ApAvEvTimes(3);
    ApAvPavlovEvs=(findEvents(outcomes('ForcedRewardApAv'),[216])+findEvents(outcomes('ForcedPunishApAv'),[212]))/2;
    
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


    