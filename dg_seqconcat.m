function seq = dg_seqconcat(delimiter, varargin)
%seq = dg_seqconcat(delimiter, varargin)
% Each argument in <varargin> is a cell vector of column vectors specifying
% sequences of symbols.  Like Perl's "join", concatenates all the
% sequences with the value <delimiter> placed in between successive
% sequence pairs.  It is an error if <delimiter> is not of a compatible
% type with the elements of the sequences.
%
% Written 9-Mar-2005 by Dan Gibson.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    error('dg_seqconcat:noseqs', ...
        'You must specify a delimiter and at least one set of sequences.' );
end

seq = varargin{1}{1};
for argnum = 1:length(varargin)
    for seqnum = 1:length(varargin{argnum})
        if seqnum > 1 || argnum > 1
            seq = [ seq; delimiter; varargin{argnum}{seqnum} ];
        end
    end
end
