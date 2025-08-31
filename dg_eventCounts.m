function varargout = dg_eventCounts(filename, varargin)
%[TTL, count, Hdr] = dg_eventCounts
%[TTL, count, Hdr] = dg_eventCounts(filename)
% If none of the outputs is used, displays the events file header and a
% table of event IDs and counts in the command window.  Assigning at least
% one output suppresses the display.
%INPUT
% <filename>:  absolute or relative pathname of file to read; if <filename>
%   is empty or not given, then a GUI is presented to choose the file.
%OUTPUTS
% <TTL>:  TTL IDs.  These are returned as signed integers, with 2^15 bit 
%	(strobe/sync) propagated to the left.  To extract the lower 15
%	bits, add 2^15 to negative values to yield positive integers.  There is
%	one element for each unique ID.
% <count>:  a vector of the same length as <TTL>, containing the number of
%   instances of each value of <TTL>.
% <Hdr>:  the events file header
%OPTIONS
% 'strobe' - reports only the counts for event IDs that have the strobe bit
%   (highest-order or sign bit) set, and extracts the value of the
%   remaining (non-strobe) bits.
% 'extract' - adds 2^15 to negative values to yield all positive integers.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin == 0
    filename = '';
end
extractflag = false;
strobeflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'extract'
            extractflag = true;
        case 'strobe'
            strobeflag = true;
        otherwise
            error('dg_eventCounts:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if isempty(filename)
    [FileName,PathName] = uigetfile({'*.nev'
        '*.dat'
        '*.*' }, ...
        'Choose events file');
    filename = fullfile(PathName, FileName);
end

[TTLdata, Hdr] = Nlx2MatEV(filename, [0, 0, 1, 0, 0], 1, 1);

if isempty(Hdr{end}) || ~isempty(regexp(Hdr{end}, '^\s*$', 'once' ))
    Hdr(end) = [];
else
    warning('dg_eventCounts:badHdr', ...
        'Header was not properly terminated with an empty line' );
end

if strobeflag
    TTLdata = TTLdata(TTLdata<0) + 2^15;
end

if extractflag
    TTLdata(TTLdata<0) = TTLdata(TTLdata<0) + 2^15;
end

TTL = unique(TTLdata);
count = zeros(size(TTL));
for TTLidx = 1:length(TTL)
    count(TTLidx) = sum(TTLdata == TTL(TTLidx));
end

if nargout == 0
    for Hdridx = 1:length(Hdr)
        disp(Hdr{Hdridx});
    end
    disp(sprintf('\n   TTLid  count'));
    for TTLidx = 1:length(TTL)
        disp(sprintf('%8.0f  %5.0f', TTL(TTLidx), count(TTLidx)));
    end
else
    varargout = {TTL, count, Hdr};
end

