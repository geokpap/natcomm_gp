function seq = dg_seqfileconcat(mydir, delimiter, varargin)
%seq = dg_seqfileconcat(mydir, delimiter,, varargin)
% Applies dg_seqconcat to all the sets of sequences in all the MAT-files
% whose names are listed in <varargin> located in directory <mydir>.  It is
% assumed that each file contains exactly one set of sequences.  The
% variable name given to the sequence set in the MAT-file is ignored.  If
% <mydir> is empty, the current working directory is used.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    error('dg_seqfileconcat:nofiles', ...
        'You must specify a directory and a list of filenames.' );
end

if isempty(mydir)
    mydir = pwd;
end
if exist(mydir) ~= 7
    error('dg_seqfileconcat:baddir', ...
        '%s is not a directory', mydir );
end

for argnum = 1:length(varargin)
    % read each sequence set from its file and save it
    S = load(fullfile(mydir, varargin{argnum}));
    fields = fieldnames(S);
    seqsets{argnum} = S.(fields{1});
    clear S;
end
% then submit whole mess to dg_seconcat
seq = dg_seqconcat(delimiter, seqsets{:});

