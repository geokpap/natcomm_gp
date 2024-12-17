load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\matrix_withDebLick.mat')
%concatZoomed=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\Lick\debbieConcatZoomed.mat');
load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\prezAnnex5Zoomed.mat')
prezA4=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\prezAnnex4BehavioralZoomed.mat');
prezA5=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\prezAnnex5BehavioralZoomed.mat');

prezA4Sessions=keys(prezA4.zoomedData);
for s=1:length(prezA4Sessions)
    zoomedData(prezA4Sessions{s})=prezA4.zoomedData(prezA4Sessions{s});
end
prezA5Sessions=keys(prezA5.zoomedData);
for s=1:length(prezA5Sessions)
    zoomedData(prezA5Sessions{s})=prezA5.zoomedData(prezA5Sessions{s});
end

%matrixData=newMatrix;

for i=132:435
    sessionName=matrixData{i,2};
    choiceMatrix=matrixData{i,3};
    pavlovMatrix=matrixData{i,4};
    try
        %{
        if contains(sessionName,'concat')
            sessionData=concatZoomed.zoomedData(sessionName);
        else
            sessionData=zoomedData(sessionName);
        end
        %}
        sessionData=zoomedData(sessionName);
    catch
        choiceMatrix(1:end,9)=NaN;
        pavlovMatrix(1:end,9)=NaN;
        matrixData{i,3}=choiceMatrix;
        matrixData{i,4}=pavlovMatrix;
        continue
    end %matrixData(151:end,:) annex4matrix(3:4,:)
    %load data
    approachData=sessionData('ApproachApAv');
    avoidData=sessionData('AvoidApAv');
    rewardData=sessionData('ForcedRewardApAv');
    airpuffData=sessionData('ForcedPunishApAv');
    %add licking to choice data
    lickData={approachData(:).Lick};
    means=mean(vertcat(lickData{:}).');
    choiceMatrix(1:length(approachData),9)=means;
    lickData={avoidData(:).Lick};
    means=mean(vertcat(lickData{:}).');
    choiceMatrix(length(approachData)+1:end,9)=means;

    %add licking to pavlovian data
    lickData={rewardData(:).Lick};
    means=mean(vertcat(lickData{:}).');
    pavlovMatrix(1:length(rewardData),9)=means;
    lickData={airpuffData(:).Lick};
    means=mean(vertcat(lickData{:}).');
    %if contains(sessionName,'2021-12-17_concat')
        %pavlovMatrix=[pavlovMatrix(1:23,:); pavlovMatrix(25:end,:)];
    %end
    pavlovMatrix(length(rewardData)+1:end,9)=means;

    matrixData{i,3}=choiceMatrix;
    matrixData{i,4}=pavlovMatrix;

end

