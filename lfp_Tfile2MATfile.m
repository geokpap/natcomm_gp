function lfp_Tfile2MATfile(sessiondir, Efilename)
%lfp_Tfile2MATfile(sessiondir, Efilename)
% Creates Multiple Cut Cluster .MAT file equivalents of all T-files from
% the session represented by <sessiondir>.  If <Efilename> is not given or
% is empty, then if 'events.nev' exists it is used; otherwise if
% 'events.dat' exists it is used; otherwise an error is raised.
% WARNING: silently overwrites any pre-existing .MAT files of the same
% name.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;
if nargin < 2 || isempty(Efilename)
    if exist(fullfile(sessiondir, 'events.nev')) == 2
        Efilename = 'events.nev';
    elseif exist(fullfile(sessiondir, 'events.dat')) == 2
        Efilename = 'events.dat';
    else
        error('lfp_Tfile2MATfile:nonev', ...
            'No events.nev or events.dat in %s', sessiondir);
    end
end

lfp_read2('preset', sessiondir, {Efilename});
[animaldir, sessionID] = fileparts(sessiondir);
Tfiles = dir(fullfile(animaldir, [sessionID '.T*']));
for tfidx = 1:length(Tfiles)
    Tfilename = Tfiles(tfidx).name;
    lfp_add('preset', animaldir, {Tfilename}, 'Rodent Clusters (*.Tnn, *.TTn)', false);
    lfp_save('preset', lfp_SpikeNames, 'mcc');
    lfp_Spikes = {};
    lfp_SpikeNames = {};
end
