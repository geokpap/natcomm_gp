% Reads the Rodent Cluster file specified by the string variable
% <filename>, and packages it into the following form. Fake absolute
% timestamps are created in the fashion of dg_makeDummyCSC. Cluster names
% of the form <animalID><sessionID>T<n>C<m> are created from the filename
% and directory name, and a spike train for each cluster is assigned to a
% variable of the same name.  Event timestamps are assigned to variables
% with names of the form Evt<n> where <n> is the event ID; the variable is
% not created until the first timestamp is encountered.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

disp('dg_RodentClust2Nex'); % this puts nothing in cmd window when called from Nex (sigh...)
[pathstr,sessionID,ext] = fileparts(filename);
[pathstr,animalID] = fileparts(pathstr);
[FileHeader, TrialData] = dg_ReadRodentFormat(filename);

% Convert event time stamps into fake-absolute timestamps and form into
% single-event-ID vectors:

eventses = cell(size(TrialData));
[eventses{:}] = deal(TrialData.events);
eventses = 100 * cell2mat(eventses');  % convert from 0.1 ms to us
eventses(:,50) = [];    % always contains garbage
% Put an ITI of 10 sec between every pair of trials:
lastevt = max(eventses, [], 2);
offset(1,1) = lastevt(1) + 5e6;
for k = 2:(numel(lastevt) - 1)
    offset(k,1) = offset(k-1,1) + lastevt(k) + 10e6;
end
offset = [0; offset];   % make life easier later

for eventID = 1:size(eventses,2)
    if any(eventses(:, eventID))
        nonzeros = find(eventses(:, eventID) ~= 0);
        command = sprintf(...
            'Evt%02d = 1e-6 * (eventses(nonzeros, eventID) + offset(nonzeros));', ...
            eventID );
        eval(command);
    end
end

% Convert spike time stamps into fake-absolute timestamps in seconds and
% assign:
for clustnum = 1:FileHeader.CSize
    digits = regexpi(ext, '[0-9]+$');
    trodenum = str2num(ext(digits:end));
    clustname = sprintf('%s%sT%dC%d', animalID, sessionID,...
        trodenum, clustnum );
    % Axe weird chars (any other than alphanumeric or underscore):
    tf = isstrprop(clustname, 'alphanum') | clustname == '_';
    clustname = clustname(tf);
    command = sprintf(...
        '%s = [];', clustname );
    eval(command);
    for trialnum = 1:length(TrialData)
        % convert from 0.1 ms to us
        spikes = 100 * (TrialData(trialnum).spikes(:,clustnum));
        nonzeros = find(spikes ~= 0);
        spikes = 1e-6 * (spikes(nonzeros) + offset(trialnum));
        command = sprintf(...
            '%s = [%s; spikes];', clustname, clustname );
        eval(command);
    end
end

clear EventTimeStamps FileHeader TrialData animalID clustname clustnum ...
    command digits eventID eventses ext k lastevt nonzeros offset ...
    pathstr sessionID spikes trialnum trialoffset trodenum
