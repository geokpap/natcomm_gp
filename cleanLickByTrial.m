load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\debbieAnnex5FullALL.mat');
load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\debbieAnnex5PlottingALL.mat')
debbieA4=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\debbieAnnex4FullALL.mat');
%prezA4=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\prezAnnex4BehavioralFullALL.mat');



%put annex 4 and annex 5 behavioral in one data structure
debbieA4Sessions=keys(debbieA4.newFullData);
for s=1:length(debbieA4Sessions)
    newFullData(debbieA4Sessions{s})=debbieA4.newFullData(debbieA4Sessions{s});
end
%{
prezA5Sessions=keys(prezBA5.newFullData);
for s=1:length(prezA5Sessions)
    newFullData(prezA5Sessions{s})=prezBA5.newFullData(prezA5Sessions{s});
end
%}
fullData=newFullData;
%find mean of 2021-2023 sessions
sessionNames=keys(fullData);
ex=fullData(sessionNames{1});
ex=ex{1};
categories={'ApproachApAv', 'AvoidApAv', 'ForcedRewardApAv', 'ForcedPunishApAv'};
means=[];
for s=1:length(sessionNames)
    sessionName=sessionNames{s};
    if contains(sessionName,'concat') | contains(sessionName,'2022-01-12_10-02-42')
        continue
    elseif contains(sessionName, '2019') | contains(sessionName, '2020')
        continue
    end
    sessionData=fullData(sessionName);
    sessionData=sessionData{1};
    if length(keys(sessionData))~=4
        continue
    end
    sessionMeans=zeros(1,4);
    for i=1:length(categories)
        catData=sessionData(categories{i});
        l=[catData.Lick];
        sessionMeans(i)=mean(l);
    end
    means=[means; sessionMeans];
end
%find standard deviation of those sessions (one value for each trial type)
STD=std(means);
meansPerCat=mean(means);
%meansPerCat=horzcat(meansPerCat,meansPerCat);
%STD=horzcat(STD,STD);
cleanedDebbieA5=cleanData(newFullData,meansPerCat,STD);

save('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\clean_debbie.mat',"cleanedDebbieA5","eventTimes","pavlovCueOff")

function cleanedAvgs=cleanData(fullData,meanLick,stdLick)
    sessions=keys(fullData);
    cleanedAvgs=containers.Map;
    outcomeLabels={'ApproachApAv', 'AvoidApAv', 'ForcedPunishApAv', 'ForcedRewardApAv'}; %{'Red','Yellow','ForcedRed','ForcedYellow'};
    for i=outcomeLabels
        cleanedAvgs(char(i))=table('Size',[length(sessions),2],'VariableTypes',["cell","cell"],'RowNames',sessions,'VariableNames',["ts","lick"]);
    end
    %rebuild same data structures only adding good trials to it
    for s=1:length(sessions)
        disp(sessions{s})
        discard=false;
        if contains(char(sessions(s)),"2019-09-17_13-40-24")
            continue
        end
        session=fullData(char(sessions(s)));
        session=session{1};
        %gather and average all the trials for one session seperated by outcome, eliminating bad ones
        %outcomes=containers.Map; 
        averages=containers.Map;
        categories={'ApproachApAv', 'AvoidApAv', 'ForcedPunishApAv', 'ForcedRewardApAv'};
        for c=1:length(categories)
            category=char(categories(c));
            categoryMean=meanLick(c);
            categorySTD=stdLick(c);
            try
                data=session(category);
            catch
                continue
            end
            categoryLickMatrix=[];
            badTrials=0;
            for i=1:length(data)
                trial=data(i);
                if mean(trial.Lick)>categoryMean+categorySTD*2 %exclude data 2std away from mean
                    badTrials=badTrials+1;
                else
                    categoryLickMatrix(end+1,:)=trial.Lick;
                    %cleanedApproach(end+1)=trial;
                    %{
                    if ismember(category,keys(outcomes))
                        old=outcomes(category);
                        old(end+1)=trial;
                        outcomes(category)=old;
                    else
                        outcomes(category)=[trial];
                    end
                    %}
                end
            end
            pctBad=badTrials/length(data);
            if pctBad>.5
                discard=true;
            end
            if ~isempty(categoryLickMatrix)
                % average the cleaned data from that category
                if length(categoryLickMatrix(:,1))==1
                    trialAvs=struct('TimeStamps',linspace(0,(length(categoryLickMatrix(1,:))-1)*.001,length(categoryLickMatrix(1,:))),'Lick',categoryLickMatrix);
                else
                    trialAvs=struct('TimeStamps',linspace(0,(length(categoryLickMatrix(1,:))-1)*.001,length(categoryLickMatrix(1,:))),'Lick',mean(categoryLickMatrix));
                end
                averages(category)=trialAvs;
            end
        end
        if discard==false
            for i=outcomeLabels
                if ismember(i,keys(averages))
                    tbl=cleanedAvgs(char(i));
                    cat=averages(char(i));
                    tbl(char(sessions(s)),:)={{cat.TimeStamps},{cat.Lick}};
                    cleanedAvgs(char(i))=tbl;
                end
            end
        end
        
        
    end
end