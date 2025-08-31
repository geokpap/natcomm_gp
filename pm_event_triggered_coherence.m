% script pm_event_triggered_coherence
% Slightly modified from Partha Mitra email
% Coherences for [-<window> <window>] sec window around event <ne> 

%$Rev: 327 $
%$Date: 2014-05-29 17:22:03 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

session = 'b10(l)acq16wSpikesv6.mat'; % session
ne=14; % event number
window = .25; % seconds before and after event
nch1=1; nch2=29; nch3=8; nch4=5; % channel numbers
ntr=40; % Number of trials

load(session);
dt=1/sampleperiod; % number of time points per second (digitisation rate)
nt=round(window*dt); % find number of samples in 1 sec interval
ts1=zeros(2*nt,ntr); ts2=ts1;ts3=ts1;ts4=ts1;
% arrays to contain time series 1 sec around the event
for j=1:ntr,
    ind=find(events{j}(:,2)==ne); % find event with label ne
    ev=round(events{j}(ind,1)*dt);
    ts1(:,j)=wavedata{j}(ev-nt:ev+nt-1,nch1);
    ts2(:,j)=wavedata{j}(ev-nt:ev+nt-1,nch2);
%     ts3(:,j)=wavedata{j}(ev-nt:ev+nt-1,nch3);
%     ts4(:,j)=wavedata{j}(ev-nt:ev+nt-1,nch4);

end
fmin=1; fmax=50;
'begin'
tic
[C12,phi,f,confC,phierr,Cerr12]=chronux_coherencyc(ts1,ts2,[2.5 3],2,dt,[fmin fmax],[2 0.05],1);
toc
'C12'
% [C23,phi,f,confC,phierr,Cerr23]=chronux_coherencyc(ts2,ts3,[2 3],2,dt,[fmin fmax],[2 0.05],1);
% toc
% 'C23'
% [C34,phi,f,confC,phierr,Cerr34]=chronux_coherencyc(ts3,ts4,[2 3],2,dt,[fmin fmax],[2 0.05],1);
% toc

% figure
% subplot(2,2,1)
% plot(f,C12); hold on; plot(f,Cerr12(1,:),'k'); plot(f,Cerr12(2,:),'k'); hold off
% axis([0 fmax 0.3 1]); 
% title(sprintf('%s coh evt %d win %d', ...
%     session, ne, window));
% ylabel(sprintf('%sx%s', filenames{nch1}, filenames{nch2}));
% subplot(2,2,2)
% plot(f,C23); hold on; plot(f,Cerr23(1,:),'k'); plot(f,Cerr23(2,:),'k'); hold off
% axis([0 fmax 0.3 1]); 
% ylabel(sprintf('%sx%s', filenames{nch2}, filenames{nch3}));
% subplot(2,2,3)
% plot(f,C34); hold on; plot(f,Cerr34(1,:),'k'); plot(f,Cerr34(2,:),'k'); hold off
% axis([0 fmax 0.3 1]); 
% ylabel(sprintf('%sx%s', filenames{nch3}, filenames{nch4}));
% subplot(2,2,4)
% plot(f,C12,'b',f,C23,'k',f,C34,'g'); 
% axis([0 fmax 0.3 1]); 

figure; 
plot(f,C12,'b'); 
% plot(f,C12,'b',f,C23,'k',f,C34,'g'); 
legend(sprintf('%sx%s', filenames{nch1}, filenames{nch2}) );
%     sprintf('%sx%s', filenames{nch2}, filenames{nch3}), ...
%     sprintf('%sx%s', filenames{nch3}, filenames{nch4}) );
axis([0 fmax 0 1]); 
title(sprintf('%s coh evt %d win %d', ...
    session, ne, window));

% fmin=1; fmax=400;
% 'begin'
% tic
% [C12,phi,f,confC,phierr,Cerr12]=chronux_coherencyc(ts1,ts2,[5 9],2,dt,[fmin fmax],[2 0.05],1);
% toc
% 'C12'
% [C23,phi,f,confC,phierr,Cerr23]=chronux_coherencyc(ts2,ts3,[5 9],2,dt,[fmin fmax],[2 0.05],1);
% toc
% 'C23'
% [C34,phi,f,confC,phierr,Cerr34]=chronux_coherencyc(ts3,ts4,[5 9],2,dt,[fmin fmax],[2 0.05],1);
% toc

% figure
% subplot(2,2,1)
% plot(f,C12); hold on; plot(f,Cerr12(1,:),'k'); plot(f,Cerr12(2,:),'k'); hold off
% ylabel(sprintf('%sx%s', filenames{nch1}, filenames{nch2}));
% axis([0 fmax 0.3 1]); 
% title(sprintf('%s coh evt %d win %d', ...
%     session, ne, window));
% subplot(2,2,2)
% plot(f,C23); hold on; plot(f,Cerr23(1,:),'k'); plot(f,Cerr23(2,:),'k'); hold off
% ylabel(sprintf('%sx%s', filenames{nch2}, filenames{nch3}));
% axis([0 fmax 0.3 1]); 
% subplot(2,2,3)
% plot(f,C34); hold on; plot(f,Cerr34(1,:),'k'); plot(f,Cerr34(2,:),'k'); hold off
% ylabel(sprintf('%sx%s', filenames{nch3}, filenames{nch4}));
% axis([0 fmax 0.3 1]); 
% subplot(2,2,4)
% plot(f,C12,'b',f,C23,'k',f,C34,'g'); 
% axis([0 fmax 0.3 1]); 
% 
% figure; 
% plot(f,C12,'b',f,C23,'k',f,C34,'g'); 
% legend(sprintf('%sx%s', filenames{nch1}, filenames{nch2}), ...
%     sprintf('%sx%s', filenames{nch2}, filenames{nch3}), ...
%     sprintf('%sx%s', filenames{nch3}, filenames{nch4}) );
% axis([0 fmax 0.3 1]); 
% title(sprintf('%s coh evt %d win %d', ...
%     session, ne, window));
% 
