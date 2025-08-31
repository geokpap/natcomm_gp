function timerange = lfp_findTimeRange(filename)
%timerange = lfp_findTimeRange(filename)
%timerange = lfp_findTimeRange
%   Gets time range in raw Neuralynx ticks spanned by CSC file <filename>.
%   If <filename> is empty or not given, then the first file that matches
%   CSCFileRegexp, CSCFileExt from lfp_getEvtIDs_* is used.  If there is no
%   such file, returns [].

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

lfp_declareGlobals;
if nargin < 1
    filename = '';
end
lfp_getEvtIDs;
if ~isempty(filename)
    CSCfilelist = {filename};
else
    DataFiles = dir(lfp_DataDir);
    CSCfilelist = {};
    for entrynum = 1:length(DataFiles)
        file = DataFiles(entrynum);
        filenameCSC = file.name;
        filenameCSC = upper(filenameCSC);
        % This is where we decide which file to open:
        if regexpi(filenameCSC, [CSCFileRegexp '\.' CSCFileExt])
            CSCfilelist = {filenameCSC};
            break
        end
    end
end
if isempty(CSCfilelist)
    timerange = [];
    return
end
pathname = fullfile(lfp_DataDir, CSCfilelist{1});
[pathstr,name,ext] = fileparts(pathname); %#ok<ASGLU>
switch (upper(ext))
    case '.MAT'
        load(pathname, '-mat');
        CSCTimeStamps = dg_Nlx2Mat_Timestamps;
        clear dg_Nlx2Mat_Timestamps;
        clear dg_Nlx2Mat_Samples;
    otherwise
        CSCTimeStamps = Nlx2MatCSC(pathname, [1, 0, 0, 0, 0], 0, 1);
end
framedur = CSCTimeStamps(2) - CSCTimeStamps(1);
minTS = CSCTimeStamps(1);
maxTS = CSCTimeStamps(end) + framedur;
timerange = [minTS maxTS];
