load('X:\georgios\UROP\Victoria\Matlab\Data\ApAp_det_rev\Megatask\debbieMatrix_zero.mat')
%load('X:\georgios\UROP\Victoria\Matlab\Data\ApAp_det_rev\Zero\HR\prezMatrix.mat')
%pSessions=keys(prezSessions);
pSessions=[];
sessions=keys(debbieSessions);
matrixData=cell(length(sessions)+length(pSessions),4);
for i=1:length(sessions)
    session=char(sessions(i));
    sessionMatrices=debbieSessions(session);
    sessionMatrix=sessionMatrices{1};
    pavlovSessionMatrix=sessionMatrices{2};
    matrixData{i,1}='Debbie';
    matrixData{i,2}=session;
    matrixData{i,3}=sessionMatrix;
    matrixData{i,4}=pavlovSessionMatrix;
end

for i=1:length(pSessions)
    session=char(pSessions(i));
    sessionMatrices=prezSessions(session);
    sessionMatrix=sessionMatrices{1};
    pavlovSessionMatrix=sessionMatrices{2};
    idx=length(sessions)+i;
    matrixData{idx,1}='Prez';
    matrixData{idx,2}=session;
    matrixData{idx,3}=sessionMatrix;
    matrixData{idx,4}=pavlovSessionMatrix;
end

columnLabels=[["Red Size","Yellow Size","1/0 for Approach/Avoidance", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Block","TrialNum"];
    ["Reward/Airpuff Size","1/0 for Forced Red/Forced Yellow","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV","Block","TrialNum"]];

save('X:\georgios\UROP\Victoria\Matlab\Data\ApAp_det_rev\Megatask\ApAv_det_rev_airpuff_zero.mat',"matrixData","columnLabels")

