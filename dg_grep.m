function dg_grep(regexp, varargin)

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

thingflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'thing'
            thingflag = true;
        otherwise
            error('funcname:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
