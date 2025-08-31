% Slightly modified from Partha Mitra email
% Spectra for [-<window> <window>] sec window around event <ne>

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

session = 'b08acq06_T6C1'; % session
ne=14; % event number
window = 1; % seconds before and after event
nch=7; % channel number
ntr=40; % Number of trials

load(session);
dt=1/sampleperiod; % number of time points per second (digitisation rate)
nt=round(window*dt); % find number of samples in 0.5 sec interval
ts=zeros(2*nt,ntr); % array to contain time series 1 sec around the event
npad=8192; 
spec=zeros(npad/2,ntr); specd=spec;
for j=1:ntr,
    ind=find(events{j}(:,2)==ne); % find event with label ne
    ev=round(events{j}(ind,1)*dt);
    ts(:,j)=wavedata{j}(ev-2*nt:ev-1,nch);
end
tsm=mean(ts,2);
tsd=ts-repmat(tsm,1,ntr);
for j=1:ntr, 
    [sp fr]=pmtm(tsd(:,j),2,npad,dt); 
    [spd fr]=pmtm(diff(tsd(:,j)),2,npad,dt); 
    spec(:,j)=sp(1:npad/2); specd(:,j)=spd(1:npad/2); 
end
fr=fr(1:npad/2);
figure;
imagesc([1:40],fr(1:1000),10*log10(spec(1:1000,:)))
axis xy
title(sprintf('%s Sp1 %s evt %d win %d', ...
    session, filenames{nch}, ne, window));

figure;
plot(fr,10*log10(mean(spec,2))); 
title(sprintf('%s Sp2 %s evt %d win %d', ...
    session, filenames{nch}, ne, window));

figure;
plot(fr(1:1000),10*log10(mean(spec(1:1000,:),2)));
title(sprintf('%s Sp3 %s evt %d win %d', ...
    session, filenames{nch}, ne, window));

figure;
imagesc([1:40],fr(35:1000),10*log10(specd(35:1000,:)))
axis xy 
title(sprintf('%s Sp4 %s evt %d win %d', ...
    session, filenames{nch}, ne, window));

figure;
plot(fr(35:1000),10*log10(mean(specd(35:1000,:),2))); 
title(sprintf('%s Sp5 %s evt %d win %d', ...
    session, filenames{nch}, ne, window));

