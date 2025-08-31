function lfp_editRodentTT(infilename, evtfilename, outfilename, varargin)
%lfp_editRodentTT(infilename, evtfilename, outfilename)
%	Reads the Neuralynx tetrode data file <infilename>, copies data to
%	Neuralynx tetrode data file <outfilename> only for the time periods
%	starting 2 s before the lfp_NominalTrialStart event and ending at
%	lfp_NominalTrialEnd.  <evtfilename> must be an events file that can be
%	read by lfp_read2, in the same directory as <infilename>.
%   <infilename> must include the full path to the session data directory;
%   <evtfilename> must not. If <infilename> had the extension '.dat', then
%   it is copied to the current working directory with the extension
%   changed to '.ntt' before reading, and the copy is deleted after.  If
%   there is already such a file in your current working directory, then
%   the copy is not done but instead a warning is issued, the file is still
%   read, and it is not deleted.
% OPTIONS
% 'skipbadrecs' - deletes records whose timestamps are out of temporal
%   order. 
% 'addperiods', evtIDs, auxdata - adds time periods defined by <evtIDs>
%   and/or <auxdata>.  <evtIDs> is a two-element vector containing a start
%   event ID and an end event ID that define each additional period to
%   include.  If <auxdata> is empty, then it is assumed that <evtfilename>
%   already contains events with the IDs given in <evtIDs>.  Otherwise,
%   <auxdata> should be a cell array whose first element is either empty or
%   is a string to be passed to lfp_save as the name for a new version of
%   the events data.  If it's empty, no file is saved; otherwise, the
%   '.evtsav' file extension will be automatically added. The second
%   element in the cell array should be a two-column list of timestamps
%   whose first column will be used to insert <evtIDs(1)> and whose second
%   column will be used to insert <evtIDs(2)>.  It is assumed that the time
%   periods for trials and additional periods are strictly disjoint, but if
%   this is not the case then a warning is issued.

%$Rev: 333 $
%$Date: 2014-10-16 18:56:27 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

argnum = 1;
auxdata = [];
evtIDs = [];
skipbadrecsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'addperiods'
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('lfp_editRodentTT:addperiods1', ...
                    '''addperiods'' requires additional args <evtIDs> and <auxdata>');
            end
            evtIDs = varargin{argnum};
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('lfp_editRodentTT:addperiods2', ...
                    '''addperiods'' requires additional args <evtIDs> and <auxdata>');
            end
            auxdata = varargin{argnum};
            skipbadrecsflag = true;
        case 'skipbadrecs'
            skipbadrecsflag = true;
        otherwise
            error('lfp_editRodentTT:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

lfp_declareGlobals;
lfp_getEvtIDs;
disp(['Using ' lfp_SetupName ' setup']);
[pathstr, name, ext] = fileparts(infilename);
deltmpfile = false;
if strcmpi(ext, '.dat')
    tmpfilename = [name '.ntt'];
    if exist(tmpfilename, 'file')
        warning('lfp_editRodentTT:gotfile', ...
            'The file %s already exists in your current working directory; skipping copy', ...
            tmpfilename );
    else
        copyfile(infilename, tmpfilename);
        deltmpfile = true;
    end
    infilename = tmpfilename;
end
lfp_read2('preset', pathstr, {evtfilename, ''});
% selectedtime contains start times in col 1, end times in col 2:
selectedtime(:, 1) = (lfp_Events(lfp_TrialIndex(:,1), 1) - 2) * 1e6;
selectedtime(:, 2) = (lfp_Events(lfp_TrialIndex(:,2), 1)) * 1e6;
if ~isempty(evtIDs)
    if ~isempty(auxdata)
        if ~isempty(auxdata{1}) && ~ischar(auxdata{1})
            error('lfp_editRodentTT:badname', ...
                'The <auxdata> arg to ''addperiods'' must contain a string as its first element');
        end
        if size(auxdata{2},2) ~= 2 || size(auxdata{2},1) == 0
            error('lfp_editRodentTT:badTS', ...
                'The <auxdata> arg to ''addperiods'' must contain a two-column array as its second element');
        end
        lfp_createEvents(@(a) {a(:,1) a(:,2)}, evtIDs, auxdata{2});
        if ~isempty(auxdata{1})
            lfp_save('preset', auxdata{1}, 'evt');
        end
    end
    addtime(:, 1) = (lfp_Events(lfp_Events(:,2)==evtIDs(1), 1)) * 1e6;
    if sum(lfp_Events(:,2)==evtIDs(2)) ~= size(addtime,1)
        error('lfp_editRodentTT:badevts', ...
            'There are unequal numbers of start and end events for ''addperiods''');
    end
    addtime(:, 2) = (lfp_Events(lfp_Events(:,2)==evtIDs(2), 1)) * 1e6;
    for k = 1:size(addtime,1)
        msg = '';
        if any( selectedtime(:,1) <= addtime(k,1) ...
                & selectedtime(:,2) >= addtime(k,1) )
            % start of addtime is inside selectedtime; delete from addtime
            % and set endpoint in selectedtime to the later of the two
            % endpoints
            msg = sprintf( ...
                'addperiods start time %.6f is inside a trial', ...
                addtime(k,1)*1e-6 );
        end
        if any( selectedtime(:,1) <= addtime(k,2) ...
                & selectedtime(:,2) >= addtime(k,2) )
            % end of addtime is inside selectedtime; delete from addtime
            % and set endpoint in selectedtime to the later of the two
            % endpoints
            msg = sprintf( ...
                'addperiods end time %.6f is inside a trial', ...
                addtime(k,2)*1e-6 );
        end
        if ~isempty(msg)
            warning('lfp_editRodentTT:overlap1', '%s', msg);
        end
    end
    for k = 1:size(selectedtime,1)
        msg = '';
        if any( addtime(:,1) <= selectedtime(k,1) ...
                & addtime(:,2) >= selectedtime(k,1) )
            % start of selectedtime is inside addtime; delete from selectedtime
            % and set endpoint in addtime to the later of the two
            % endpoints
            msg = sprintf( ...
                'trial start time %.6f is inside an additional period', ...
                selectedtime(k,1)*1e-6 );
        end
        if any( addtime(:,1) <= selectedtime(k,2) ...
                & addtime(:,2) >= selectedtime(k,2) )
            % end of selectedtime is inside addtime; delete from selectedtime
            % and set endpoint in addtime to the later of the two
            % endpoints
            msg = sprintf( ...
                'trial end time %.6f is inside an additional period', ...
                selectedtime(k,2)*1e-6 );
        end
        if ~isempty(msg)
            warning('lfp_editRodentTT:overlap2', '%s', msg);
        end
    end
    selectedtime = [selectedtime; addtime];
    selectedtime = sortrows(selectedtime);
end

% process file in batches no bigger than reclimit
reclimit = 200000;
[allTS, NlxHeader] = Nlx2MatSpike_411(infilename, [1,0,0,0,0], 1, 1);
if any(allTS(2:end) < allTS(1:end-1))
    badrecs = find(allTS(2:end) < allTS(1:end-1));
    msgfmt = '%.0f disordered spike time(s) starting at record %d, TS %.6f:\n%s';
    if skipbadrecsflag
        numTS = numel(allTS);
        warning('lfp_editRodentTT:badfile', ...
            [msgfmt, '\nDeleting bad records %s.'], ...
            length(badrecs), ...
            badrecs(1)+1, allTS(badrecs(1)+1)*1e-6, infilename, ...
            dg_canonicalSeries(badrecs(1)+1) );
        if length(badrecs) > 5
            warning('lfp_editRodentTT:worsefile', ...
                'ATTENTION: this deletion involves %.0f spike records!', ...
                length(badrecs) );
        end
        while ~isempty(badrecs)
            allTS(badrecs+1) = [];
            badrecs = find(allTS(2:end) < allTS(1:end-1));
        end
        lfp_log(sprintf('Deleted %d records from %s.', ...
            numTS - numel(allTS), infilename));
    else
        evtmsg = sprintf('lfp_Events:\n%14s\tevtID', 'TimeStamp');
        for k = 1:5
            evtmsg = sprintf('%s\n%14.6f\t%5.0f', ...
                evtmsg, lfp_Events(k,:));
        end
        badTSmsg = 'Bad Spike Times:';
        for k = 1:length(badrecs)
            badTSmsg = sprintf('%s\n%14.6f', ...
                badTSmsg, allTS(badrecs(k)+1)*1e-6);
        end
        error('lfp_editRodentTT:badfile2', ...
            msgfmt, length(badrecs), ...
            badrecs(1)+1, allTS(badrecs(1)+1)*1e-6, ...
            sprintf('%s\n%s\n%s', infilename, evtmsg, badTSmsg) );
    end
end
numrecs = length(allTS);
numbatches = ceil(numrecs/reclimit);
for batchnum = 1:numbatches
    firstrec = (batchnum-1) * reclimit + 1;
    lastrec = min(numrecs, batchnum * reclimit);
    [ScNumbers, CellNumbers, Params, DataPoints] = ...
        Nlx2MatSpike_411(infilename, [0,1,1,1,1], 0, ...
        4, [allTS(firstrec) allTS(lastrec)]);
    batchTS = allTS(firstrec:lastrec);

    % This could be made faster by using the fact that allTS is sorted:
    selectedspikes = false(size(batchTS));
    for trial = 1:size(selectedtime,1)
        selectedspikes( ...
            batchTS >= selectedtime(trial, 1) ...
            & batchTS <= selectedtime(trial, 2) ) = true;
    end

    numspikes = sum(selectedspikes);
    if numspikes
        if ~ispc
            error('lfp_editRodentTT:notPC', ...
                'Neuralynx format files can only be created on Windows machines.');
        end
        Mat2NlxTT(outfilename, batchnum>1, 1, 1, numspikes, ...
            [1 1 1 1 1 1], ...
            batchTS(selectedspikes), ScNumbers(selectedspikes), ...
            CellNumbers(selectedspikes), Params(:,selectedspikes), ...
            DataPoints(:,:,selectedspikes), NlxHeader );
    end
end

if deltmpfile
    delete(tmpfilename);
end

    
