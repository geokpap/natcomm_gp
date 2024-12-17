%create empty map
prezSessions=containers.Map;
%load outcomes
load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\prezAnnex4Full.mat')
%load('X:\georgios\UROP\Victoria\Matlab\Data\ApAvApAp\prezA4Rewards.mat')
%load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\rxnTimes\prezA4.mat')
%iterate through each session
sessions=keys(fullData);
for s=1:length(sessions)
    sessionName=char(sessions(s));
    outcomes=fullData(sessionName);
    approachData=outcomes('ApproachApAv');
    avoidData=outcomes('AvoidApAv');
    sessionMatrix=zeros(length(approachData)+length(avoidData),8);%create empty matrix of size [numapproach+num avoid trials, 7]
    %approachTrials
    appRxnTimes=reactionTimes(sessionName,'ApproachApAv');
    appRxnTimes=cell2mat(appRxnTimes.ApproachApAv);
    for i=1:length(approachData)%iterate through each trial
        trial=approachData(i);
        if length(trial.Rewards)==1 %if there is only one number it was the same for approach and avoid
            num=trial.Rewards;
            trial.Rewards=[num num];
        elseif length(trial.Rewards)~=2
            continue
        end
        sessionMatrix(i,1)=trial.Rewards(1);
        sessionMatrix(i,2)=trial.Rewards(2);
        sessionMatrix(i,3)=1;
        sessionMatrix(i,4)=mean(trial.Heartrate);
        sessionMatrix(i,5)=mean(zscore(trial.Heartrate));
        sessionMatrix(i,6)=mean(trial.HRV);
        sessionMatrix(i,7)=mean(zscore(trial.HRV));
        sessionMatrix(i,8)=appRxnTimes(i);
    end
    %avoid trials
    avRxnTimes=reactionTimes(sessionName,'AvoidanceApAv');
    avRxnTimes=cell2mat(avRxnTimes.AvoidanceApAv);
    for i=1:length(avoidData)%iterate through each trial
        trial=avoidData(i);
        if length(trial.Rewards)==1 %if there is only one number it was the same for approach and avoid
            num=trial.Rewards;
            trial.Rewards=[num num];
        elseif length(trial.Rewards)~=2
            continue
        end
        idx=length(approachData)+i;
        sessionMatrix(idx,1)=trial.Rewards(1);% save reward,punish,outcome,mean hr, mean normalized hr,meanhrv, mean normalized hrv
        sessionMatrix(idx,2)=trial.Rewards(2);
        sessionMatrix(idx,3)=0;
        sessionMatrix(idx,4)=mean(trial.Heartrate);
        sessionMatrix(idx,5)=mean(zscore(trial.Heartrate));
        sessionMatrix(idx,6)=mean(trial.HRV);
        sessionMatrix(idx,7)=mean(zscore(trial.HRV));
        sessionMatrix(idx,8)=avRxnTimes(i);
    end
    rewardData=outcomes('ForcedRewardApAv');
    punishData=outcomes('ForcedPunishApAv');
    pavlovSessionMatrix=zeros(length(rewardData)+length(punishData),8);
    %forced reward Trials
    for i=1:length(rewardData)%iterate through each trial
        trial=rewardData(i);
        if length(trial.Rewards)~=1
            continue
        end
        pavlovSessionMatrix(i,1)=trial.Rewards;
        pavlovSessionMatrix(i,2)=1;
        pavlovSessionMatrix(i,3)=mean(trial.Heartrate);
        pavlovSessionMatrix(i,4)=median(trial.Heartrate);
        pavlovSessionMatrix(i,5)=mean(zscore(trial.Heartrate));
        pavlovSessionMatrix(i,6)=mean(trial.HRV);
        pavlovSessionMatrix(i,7)=median(trial.HRV);
        pavlovSessionMatrix(i,8)=mean(zscore(trial.HRV));
    end
    %forced airpuff trials
    for i=1:length(punishData)%iterate through each trial
        trial=punishData(i);
        if length(trial.Rewards)~=1
            continue
        end
        idx=length(rewardData)+i;
        pavlovSessionMatrix(idx,1)=trial.Rewards;
        pavlovSessionMatrix(idx,2)=0;
        pavlovSessionMatrix(idx,3)=mean(trial.Heartrate);
        pavlovSessionMatrix(idx,4)=median(trial.Heartrate);
        pavlovSessionMatrix(idx,5)=mean(zscore(trial.Heartrate));
        pavlovSessionMatrix(idx,6)=mean(trial.HRV);
        pavlovSessionMatrix(idx,7)=median(trial.HRV);
        pavlovSessionMatrix(idx,8)=mean(zscore(trial.HRV));
    end
    prezSessions(char(sessions(s)))={sessionMatrix,pavlovSessionMatrix};
end
%debbieSessions=prezSessions;
save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\debbieA5Matrix.mat',"prezSessions")
