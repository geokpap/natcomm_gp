function [TS, x, y, Hdr] = dg_readVT(filename, varargin)
% Does not read return values that are not used.  Note that this is based
% purely on the number of return values, so values whose return is
% suppressed by using '~' *do* still get read.
% OPTIONS
%   'header' - reads only header, returns empty for TS and x, y
%   'matlab' - uses dg_Nlx2MatVT instead of Neuralynx function.  Note that
%       dg_Nlx2MatVT does not offer the Neuralynx "Extraction Mode"
%       or header-only features.
%   'mode', modenum, modearg - invokes the Neuralynx "Extraction Mode"
%       specified by <modenum> (default is 1).  In keeping with the new
%       convention of Nlx library v5.0.1, the first record is record 1,
%       whereas in releases through v4.1.3 it was record 0.
%         1. Extract All
%         2. Extract Record Index Range
%         3. Extract Record Index List
%         4. Extract Timestamp Range
%         5. Extract Timestamp List

%$Rev: 190 $
%$Date: 2014-02-06 20:37:55 -0500 (Thu, 06 Feb 2014) $
%$Author: dgibson $

TS = [];
x = [];
y = [];
selectary = [1, 0, 0, 0, 0, 0];
if nargout >= 2
    selectary(2) = 1;
end
if nargout >= 3
    selectary(3) = 1;
end
if nargout >= 4
    readheader = 1;
else
    readheader = 0;
end
headeronly = false;
matlabflag = false;
modenum = 1;
modearg = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'header'
            headeronly = true;
        case 'matlab'
            matlabflag = true;
        case 'mode'
            argnum = argnum + 1;
            modenum = varargin{argnum};
            argnum = argnum + 1;
            modearg = varargin{argnum};
        otherwise
            error('dg_readCSC:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if matlabflag
    [ts x y a targ pts hdr] = dg_Nlx2MatVT(filename); %#ok<ASGLU>
    TS = double(ts);
    x = double(x);
    y = double(y);
    hdrbt = uint8(zeros(2*length(hdr),1));
    hdrbt(1:2:end) = bitand(255, hdr);
    hdrbt(2:2:end) = bitshift(hdr, -8);
    Hdr = {''};
    for k = 1:length(hdrbt)
        if ~ismember(hdrbt(k), [13, 10, 0])
            Hdr{end}(1,end+1) = char(hdrbt(k));
        else
            % found <CR>, <LF>, or <NUL>
            if ~isempty(Hdr{end})
                Hdr(end+1,1) = {''};
            end
        end
    end
else
    [p, n, ext] = fileparts(filename); %#ok<ASGLU>
    if strcmpi(ext, '.dat')
        warning('dg_readCSC:ext', ...
            'This Nlx function requires .nvt extension; making temporary copy of %s.', ...
            filename );
        tempfn = [tempname '.nvt'];
        while exist(tempfn) %#ok<EXIST>
            tempfn = [tempname '.nvt'];
        end
        dg_copyfile(filename, tempfn);
        file2read = tempfn;
    else
        file2read = filename;
    end
    
    if ispc
        if ismember(modenum, [2 3])
            % v4.1.1 uses 0 to denote the first record
            modearg = modearg - 1;
        end
        if headeronly
            Hdr = Nlx2MatVT_411(file2read, [0, 0, 0, 0, 0], 1, 1);
        else
            if nargout > 3
                [TS, x, y, Hdr] = Nlx2MatVT_411(file2read, selectary, ...
                    readheader, modenum, modearg);
            elseif nargout == 3
                [TS, x, y] = Nlx2MatVT_411(file2read, selectary, ...
                    readheader, modenum, modearg);
            elseif nargout == 2
                [TS, x] = Nlx2MatVT_411(file2read, selectary, ...
                    readheader, modenum, modearg);
            else
                % nargout must be 1
                TS = Nlx2MatVT_411(file2read, selectary, ...
                    readheader, modenum, modearg);
            end
        end
    elseif ismac || isunix
        if ismember(modenum, [2 3])
            % The unix version uses 0 to denote the first record
            modearg = modearg - 1;
        end
        if headeronly
            Hdr = Nlx2MatVT_v3(file2read, [0, 0, 0, 0, 0], 1, 1);
        else
            if nargout > 3
                [TS, x, y, Hdr] = Nlx2MatVT_v3(file2read, selectary, ...
                    readheader, modenum, modearg);
            elseif nargout == 3
                [TS, x, y] = Nlx2MatVT_v3(file2read, selectary, ...
                    readheader, modenum, modearg);
            elseif nargout == 2
                [TS, x] = Nlx2MatVT_v3(file2read, selectary, ...
                    readheader, modenum, modearg);
            else
                % nargout must be 1
                TS = Nlx2MatVT_v3(file2read, selectary, ...
                    readheader, modenum, modearg);
            end
        end
        if strcmpi(ext, '.dat')
            dg_forcedeletefile(tempfn);
        end
    else
        error('dg_readCSC:arch', ...
            'Unrecognized computer platform');
    end
end
if exist('Hdr', 'var') && ...
        (isempty(Hdr{end}) || ~isempty(regexp(Hdr{end}, '^\s*$', 'once' )))
    Hdr(end) = [];
end

