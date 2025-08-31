function trialdata = aggroTrialPupil(rootdir)
%INPUTS
% rootdir: absolute pathname of root dir of 'trialByTrialPupilData' output
%   data tree.
%OUTPUT
% trialdata: First column is the animal name, 'Prez' or 'Debbie'.  Second
%   column is the name of the session, e.g. '2018-11-09_18-16-35fixed'.
%   Third column is the choice trials matrices, in format [size of red,
%   size of yellow, isAp, avg pupil diam,  median pupil diam, trialTS],
%   where <trialTS> denotes the timestamp in seconds of the trial's
%   lfp_NominalTrialStart event.  Fourth column is the Pavlovian trials
%   matrices, in format [size of bar, isrwd, avg pupil diam,  median pupil
%   diam, trialTS].

% Find 'trialByTrialPupilData' output files:
matfilename = 'trialByTrialPupilData.mat';
cmdstr = sprintf('find %s -name %s', rootdir, matfilename);
fprintf('Finding files...\n');
[status, cmdout] = system(cmdstr);
if status
    error('aggroTrialPupil:oops', ...
        'Failure in ''find''.');
end
matfiles = strsplit(cmdout, sprintf('\n')); %#ok<SPRINTFN>
matfiles(cellfun(@isempty, matfiles)) = [];

% Aggregate into <trialdata>, <forceddata>:
trialdata = cell(length(matfiles), 4);
for fidx = 1:length(matfiles)
    pathparts = strsplit(matfiles{fidx}, '/');
    if isequal(pathparts{end}, matfilename)
        pathparts(end) = [];
    else
        error('aggroTrialPupil:oops', ...
            'This can''t happen.');
    end
    if ~isequal(pathparts{end}(1:2), '20')
        error('aggroTrialPupil:sessionID', ...
            '%s does not end with sessionID', matfiles{fidx});
    end
    trialdata{fidx, 2} = pathparts{end};
    done = false;
    for partidx = length(pathparts) - 1 : -1 : 1
        if ismember(pathparts{partidx}, {'Prez' 'Debbie'})
            trialdata{fidx, 1} = pathparts{partidx};
            done = true;
        end
    end
    if ~done
        error('aggroTrialPupil:animalID', ...
            '%s does not contain an animal ID', matfiles{fidx});
    end
    x = load(matfiles{fidx});
    trialdata{fidx, 3} = x.choicedata;
    trialdata{fidx, 4} = x.forceddata;
end
