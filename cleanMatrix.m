load('X:\georgios\UROP\Victoria\Matlab\Data\ApAp_det_rev\Megatask\ApAv_det_rev_zero_matrix_withLick.mat')

for i=1:3
    sessionMatrix=cell2mat(matrixData(i,3));
    tf=sessionMatrix(:,4)<100 | sessionMatrix(:,4)>140;
    sessionMatrix(tf,:)=[];
    matrixData{i,3}=sessionMatrix;

    sessionMatrixPavlov=cell2mat(matrixData(i,4));
    tf=sessionMatrixPavlov(:,3)<100 | sessionMatrixPavlov(:,3)>140;
    sessionMatrixPavlov(tf,:)=[];
    matrixData{i,4}=sessionMatrixPavlov;
end
%{
columnLabels=[["Red","Yellow","1/0 for Red/Yellow", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Block","TrialNum","Lick"];
    ["Reward Size","1/0 for Forced Red/Forced Yellow","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV","Block","TrialNum","Lick"]];
%}
%{
columnLabels=[
    ["Reward Size","Airpuff Size","1/0 for Approach/Avoidance", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Lick"];
    ["Reward Size","1/0 for Forced Reward/Forced Airpuff","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV","Lick"]];

columnLabels=[["Red","Yellow","1/0 for Red/Yellow", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Lick"];
    ["Reward Size","1/0 for Forced Red/Forced Yellow","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV","Lick"];
    ["Reward Size","Airpuff Size","1/0 for Approach/Avoidance", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Lick"];
    ["Reward Size","1/0 for Forced Reward/Forced Airpuff","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV"],"Lick"];
%}
%{
columnLabels=[
    ["Red Size","Yellow Size","1/0 for Red/Yellow", "Mean HR", "Z-score HR","Mean HRV", "Z-score HRV", "Reaction Time","Lick"];
    ["Reward Size","1/0 for Forced Red/Forced Yellow","Mean HR","Median HR","Z-score HR","Mean HRV","Median HRV","Z-score HRV","Lick"]];
%}
%save('X:\georgios\UROP\Victoria\Matlab\Data\ApAvApAp\Lick\nonclean_ApAvApAp.mat',"columnLabels","matrixData")