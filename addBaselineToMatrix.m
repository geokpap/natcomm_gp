%load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\lick\prezAnnex5FullALL.mat')
%concatNewFullData=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\debbieConcatFullALL.mat');
%load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\lick\debbieAnnex5Zoomed.mat')
load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\lick\prezAnnex4FullALL.mat');
%load('X:\georgios\Analysis\Data for Heatmap Analysis\ApAv\nonclean_ApAv_withLick.mat');
sessionNames=keys(newFullData);

%TODO: run concat sessions on sort for plotting and add them to script 

sessionBaselines=containers.Map;
for s=1:length(sessionNames)
    %sessionDataZoomed=zoomedData(sessionNames{s});
    sessionDataFull=newFullData(sessionNames{s});
    fullOutcomes=sessionDataFull{1};
    %find baseline from full data
    approachData=fullOutcomes('ApproachApAv');
    avoidData=fullOutcomes('AvoidApAv');
    baselines={0 0 0 0};
    approachBaselines=zeros(length(approachData),1);
    for i=1:length(approachData)
        trial=approachData(i);
        %get 2000-4000 for baseline and store in matrix
        approachBaselines(i)=mean(trial.Lick(2000:4000));
    end
    baselines{1}=approachBaselines;
    %get baseline for all trial types
    avoidBaselines=zeros(length(avoidData),1);
    for i=1:length(avoidData)
        trial=avoidData(i);
        avoidBaselines(i)=mean(trial.Lick(2000:4000));
    end
    baselines{2}=avoidBaselines;
    %pavlovian reward
    rewardData=fullOutcomes('ForcedRewardApAv');
    airpuffData=fullOutcomes('ForcedPunishApAv');
    rewardBaselines=zeros(length(rewardData),1);
    for i=1:length(rewardData)
        trial=rewardData(i);
        rewardBaselines(i)=mean(trial.Lick(2000:4000));
    end
    baselines{3}=rewardBaselines;
    %pavlovian airpuff
    airpuffBaselines=zeros(length(airpuffData),1);
    for i=1:length(airpuffData)
        trial=airpuffData(i);
        airpuffBaselines(i)=mean(trial.Lick(2000:4000));
    end
    baselines{4}=airpuffBaselines;
    sessionBaselines(sessionNames{s})=baselines;
end

%create new column in matrix that is col 9-baseline
for i=210:211
    sessionName=matrixData{i,2};
    choiceMatrix=matrixData{i,3};
    pavlovMatrix=matrixData{i,4};
    try
        baseline=sessionBaselines(sessionName);
        sessionData=newFullData(sessionName);
        sessionData=sessionData{1};
    catch
        choiceMatrix(1:end,10)=NaN;
        pavlovMatrix(1:end,10)=NaN;
        matrixData{i,3}=choiceMatrix;
        matrixData{i,4}=pavlovMatrix;
        continue
    end %matrixData(151:end,:) annex4matrix(3:4,:)
    %load data
    approachData=sessionData('ApproachApAv');
    avoidData=sessionData('AvoidApAv');
    rewardData=sessionData('ForcedRewardApAv');
    airpuffData=sessionData('ForcedPunishApAv');
    %subtract baseline from licking data
    choiceMatrix(1:length(approachData),10)=choiceMatrix(1:length(approachData),9)-baseline{1};
    choiceMatrix(length(approachData)+1:end,10)=choiceMatrix(length(approachData)+1:end,9)-baseline{2};

    %add licking to pavlovian data
    pavlovMatrix(1:length(rewardData),10)=pavlovMatrix(1:length(rewardData),9)-baseline{3};
    pavlovMatrix(length(rewardData)+1:end,10)=pavlovMatrix(length(rewardData)+1:end,9)-baseline{4};

    matrixData{i,3}=choiceMatrix;
    matrixData{i,4}=pavlovMatrix;

end



