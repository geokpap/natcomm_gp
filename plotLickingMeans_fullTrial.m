%debA4=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAvApAp\Lick\clean_debbieAnnex4Plotting3.mat');
debA5=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\clean_debbie.mat');
%debConcat=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv\lick\clean_debbieConcatPlotting2.mat');
%prezA4=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAvApAp\lick\clean_prezAnnex4Plotting3.mat');
prezA5=load('X:\georgios\UROP\Victoria\Matlab\Data\ApAv_rev_only\Lick\clean_prez.mat');


%{
for i=1:length(categories)
    debA4Data=debA4.avgTables(char(categories(i)));
    debA4.avgTables(char(categories(i)))=rmBadSessions(debA4Data,debA4Data.Properties.RowNames);
    debA5Data=debA5.avgTables(char(categories(i)));
    debA5.avgTables(char(categories(i)))=rmBadSessions(debA5Data,debA5Data.Properties.RowNames);
    prezA4Data=prezA4.avgTables(char(categories(i)));
    prezA4.avgTables(char(categories(i)))=rmBadSessions(prezA4Data,prezA4Data.Properties.RowNames);
    prezA5Data=prezA5.avgTables(char(categories(i)));
    prezA5.avgTables(char(categories(i)))=rmBadSessions(prezA5Data,prezA5Data.Properties.RowNames);
end
%}
%combinedDebbie=concatData(debA4.cleanedDebbieA4,debA5.cleanedDebbieA5);
%combinedPrez=concatData(prezA4.cleanedPrezA4,prezA5.cleanedPrezA5);
combinedDebbie=debA5.cleanedDebbieA5;
combinedPrez=prezA5.cleanedPrezA5;
eventTimes=mean([debA5.eventTimes;prezA5.eventTimes]);
pavlovCueOff=mean([debA5.pavlovCueOff;prezA5.pavlovCueOff]);
%{
cats=keys(debConcat.cleanedDebbieConcat);
for c=1:length(cats)
    catData=debConcat.cleanedDebbieConcat(cats{c});
    allData=combinedDebbie(cats{c});
    sesns=catData.Properties.RowNames;
    for i=1:length(sesns)
        allData(sesns{i},:)=catData(sesns{i},:);
    end
end
%}
avgsOfSessions=containers.Map;
categories=["ApproachApAv" "AvoidApAv" "ForcedRewardApAv" "ForcedPunishApAv"]; %["Red","Yellow","ForcedRed","ForcedYellow"];
for o=categories
    outcome=combinedPrez(char(o));
    outcomeDebbie=combinedDebbie(char(o));
    %get avg length of trial using ts and avg # of trials using hrv over all the sessions
    %avg all the sessions and create new time vector
    tLen=0;
    numSessions=length(outcome.ts);
    nonempty=0;
    for i=1:numSessions
        time=cell2mat(outcome.ts(i));
        tLen=tLen+length(time);
        if ~isempty(time)
            nonempty=nonempty+1;
        end
    end
    numSessionsDeb=length(outcomeDebbie.ts);
    nonemptyDeb=0;
    tLenDeb=0;
    for i=1:numSessionsDeb
        time=cell2mat(outcomeDebbie.ts(i));
        tLenDeb=tLenDeb+length(time);
        if ~isempty(time)
            nonemptyDeb=nonemptyDeb+1;
        end
    end
    avTrialLen=round(((tLenDeb/nonemptyDeb)+(tLen/nonempty))/2);
    lickMatrixPrez=zeros(nonempty,avTrialLen);%create empty matrix where each row is a trial and each col is timepoint
    lickMatrixDebbie=zeros(nonemptyDeb,avTrialLen);%create empty matrix where each row is a trial and each col is timepoint
    idx=1;
    %interpolate each session to be the same length
    for i=1:numSessions
        hr=cell2mat(outcome.lick(i)); 
        if ~isempty(hr)
            lickMatrixPrez(idx,:)=interp1(hr,linspace(1,length(hr),avTrialLen)); %interpolate data to average length of each session
            idx=idx+1;
        end
    end
    %hr=cell2mat(outcome.lick(4));
    %lickMatrixPrez=interp1(hr,linspace(1,length(hr),avTrialLen)); %interpolate data to average length of each session
    idx=1;
    for i=1:numSessionsDeb
        hr=cell2mat(outcomeDebbie.lick(i));
        if ~isempty(hr)
            lickMatrixDebbie(idx,:)=interp1(hr,linspace(1,length(hr),avTrialLen)); %interpolate data to average length of each session
            idx=idx+1;
        end
    end
    
    %hrvvbefore=length(hrVector);
    %hrbefore=length(hrMatrix);
    %hrvBefore=length(hrvMatrix);

    %hrVector2=rmoutliers(hrVector,'median','ThresholdFactor',2);
    %try to clean mean and hrv data
    %hrMatrix=rmoutliers(hrMatrix);
    %hrvMatrix=rmoutliers(hrvMatrix);
    
    %{
    means=mean(hrMatrix.');
    deleted=0;
    for i=1:length(means)
        m=means(i);
        if m<90
            hrMatrix(i-deleted,:)=[];
            deleted=deleted+1;
        end
    end
    %}
    %disp('% outliers')
    %msg=['HRVs: ',int2str((hrvvbefore-length(hrVector))/hrvvbefore),' mean HR: ',int2str((hrbefore-length(hrMatrix))/hrbefore),' median HRV: ',int2str((hrvbefore-length(hrvMatrix))/hrvbefore)];
    %disp(msg)

    %prezHR=mean(lickMatrixPrez);
    prezHR=lickMatrixPrez;
    debbieHR=mean(lickMatrixDebbie);

    combinedHR=vertcat(prezHR,debbieHR);


    hr2STDd=std(lickMatrixDebbie)/sqrt(length(lickMatrixDebbie)); %calculate standard error at each point
    hr2STDp=std(lickMatrixPrez)/sqrt(length(lickMatrixPrez)); %calculate standard error at each point
    %hrv2STDp=std(hrvMatrixPrez)/sqrt(length(hrvMatrixPrez));

    %hr2STD=mean(vertcat(hr2STDd,hr2STDp));
    %hrv2STD=mean(vertcat(hrv2STDd,hrv2STDp));

    %average over all trials for that category
    %hrvout=isoutlier(hrvMatrix);
    %cleanedHRVMatrix=hrvMatrix(hrvout);
    avs=struct('TimeStamps',linspace(0,(avTrialLen-1)*.001,avTrialLen),'Lick',mean(combinedHR),'STD',hr2STDd);
    avgsOfSessions(char(o))=avs;
end

figure
hold on
d1=avgsOfSessions('ApproachApAv');
d2=avgsOfSessions('AvoidApAv');
%plot(d1.TimeStamps,d1.Lick,'color',"g")
%plot(d2.TimeStamps,d2.Lick,'color',"r")
shadedErrorBar(d1.TimeStamps,d1.Lick*1e+6,d1.STD,'lineprops','-g')
shadedErrorBar(d2.TimeStamps,d2.Lick*1e+6,d2.STD,'lineprops','-r')
xl=xline(eventTimes,'k',["Fixation Onset","Cue Onset","Peripheral Onset","Highlight On","Highlight Off","Reward/Airpuff Onset"]);
xl(1).FontSize=13;xl(2).FontSize=13;xl(3).FontSize=13;xl(4).FontSize=13;xl(5).FontSize=13;xl(6).FontSize=13;
x=[eventTimes(6) 15 15 eventTimes(6)];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
x=[0 eventTimes(1) eventTimes(1) 0];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
ylim([0 100])
xlim([.2 15])
xlabel('Time (s)')
ylabel('Licking Activity [μv]')
title('Mean Licking Activity Choice Trials')
legend('Approach','Avoidance','Location','northwest')
hold off

figure
hold on
d1=avgsOfSessions('ForcedRewardApAv');
d2=avgsOfSessions('ForcedPunishApAv');
%plot(d1.TimeStamps,d1.Lick,'color',"g")
%plot(d2.TimeStamps,d2.Lick,'color',"r")
shadedErrorBar(d1.TimeStamps,d1.Lick*1e+6,d1.STD,'lineprops','-g')
shadedErrorBar(d2.TimeStamps,d2.Lick*1e+6,d2.STD,'lineprops','-r')
xl=xline(pavlovCueOff,'k',["Fixation Onset","Cue Onset","Cue Offset","Reward/Airpuff Onset"]);
xl(1).FontSize=13;xl(2).FontSize=13;xl(3).FontSize=13;xl(4).FontSize=13;
x=[0 pavlovCueOff(1) pavlovCueOff(1) 0];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
x=[pavlovCueOff(4) 14 14 pavlovCueOff(4)];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
ylim([0 100])
xlim([.2 14])
xlabel('Time (s)')
ylabel('Licking Activity [μv]')
title('Mean Licking Activity Pavlovian Trials')
legend('Pavlovian Reward','Pavlovian Airpuff','Location','northwest')
hold off

%{
figure
hold on
d1=avgsOfSessions('Red');
d2=avgsOfSessions('Yellow');
%plot(d1.TimeStamps,d1.Lick,'color',"g")
%plot(d2.TimeStamps,d2.Lick,'color',"r")
shadedErrorBar(d1.TimeStamps,d1.Lick*1e+6,d1.STD,'lineprops','-r')
shadedErrorBar(d2.TimeStamps,d2.Lick*1e+6,d2.STD,'lineprops','-y')
xline(eventTimes,'k',["Fixation Onset","Cue Onset","Peripheral Onset","Highlight On","Highlight Off","Reward/Airpuff Onset"])
x=[eventTimes(6) 15 15 eventTimes(6)];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
x=[0 eventTimes(1) eventTimes(1) 0];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
ylim([0 100])
xlim([.2 15])
xlabel('Time (s)')
ylabel('Licking Activity [μv]')
title('Mean Licking Activity Choice Trials')
legend('Red Bar','Yellow Bar','Location','northwest')
hold off

figure
hold on
d1=avgsOfSessions('ForcedRed');
d2=avgsOfSessions('ForcedYellow');
%plot(d1.TimeStamps,d1.Lick,'color',"g")
%plot(d2.TimeStamps,d2.Lick,'color',"r")
shadedErrorBar(d1.TimeStamps,d1.Lick*1e+6,d1.STD,'lineprops','-r')
shadedErrorBar(d2.TimeStamps,d2.Lick*1e+6,d2.STD,'lineprops','-y')
xline(pavlovCueOff,'k',["Fixation Onset","Cue Onset","Cue Offset","Reward/Airpuff Onset"])
x=[0 pavlovCueOff(1) pavlovCueOff(1) 0];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
x=[pavlovCueOff(4) 14 14 pavlovCueOff(4)];
y=[0 0 100 100];
patch(x,y,'blue','FaceAlpha',.3)
ylim([0 100])
xlim([.2 14])
xlabel('Time (s)')
ylabel('Licking Activity [μv]')
title('Mean Licking Activity Pavlovian Trials')
legend('Pavlovian Red','Pavlovian Yellow','Location','northwest')
hold off
%}
%{
function combined=concatData(data1,data2)
    combined=containers.Map;
    for k=keys(data1)
        tbl1=data1(char(k));
        tbl2=data2(char(k));
        try
            tbl1('.',:)=[];
        catch
        end
        try
            tbl1('..',:)=[];
        catch
        end
        try
            tbl1('.DS_Store',:)=[];
        catch
        end
        try
            tbl1('Crashed',:)=[];
        catch
        end
        try
            tbl1('BadBehavior',:)=[];
        catch
        end
        try
            tbl1('2019-05-17_12-27-15',:)=[];
        catch
        end
        combinedTable=vertcat(tbl1,tbl2); 
        combined(char(k))=combinedTable;
    end
end

function good=rmBadSessions(data,sessionNames)
    for s=1:length(sessionNames)
        ses=char(sessionNames(s));
        if contains(ses,'2019') | contains(ses,'2020') %remove early sessions
            data(ses,:)=[];
        end
    end
    good=data;
end
%}