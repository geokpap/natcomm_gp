function dg_makeLocalAvgRefs(localavgrefs, sessiondir, prefix, destdir)
% Computes a complete set of local average reference channels for
% <sessiondir>, using dg_localAvgRef. The results are saved to files named
% <prefix><N>, where <N> = 1:length(tackle.localavgrefs).  Output files are
% written to <sessiondir>.  
% Any pre-existing file with the same names as output files are
% overwritten.
%INPUTS
% localavgrefs: a cell vector of strings or cell arrays of strings where
%   each element is a legitimate value for the <CSCfilespec> argument to
%   dg_localAvgRef that specifies the files to average together and save as
%   a "local average reference" channel.  (Note this is identical to the
%   <tackle.localavgrefs> input to lfp_fishStory.) Skips any <localavgrefs>
%   cells that are empty.
% sessiondir: absolute or relative path to directory containing CSC files.
% prefix: string to use as output filename prefix
% destdir: optional argument to specify directory for output files; default
%   value is <sessiondir>.
%OUTPUTS are only to files, named  <prefix><N>.mat, where <N> is the index
% into <localavgrefs> for the group of channels whose average in the file.
%NOTES
% dg_downsample4LocalAvgRef.m Use this when you need to downsample only the
% files that will be used in your local avg reffing
% 
% dg_localAvgRefIfNeeded.m In bulk processing, you sometimes have a mix of
% sessions that have already been local avg reffed and others that haven't.
% Use this to skip the ones that are already done.
% 
% These are the three core local avg reffing functions:
% 
% dg_localAvgRef.m Calculates one local avg ref channel from the set of CSC
% files you specify.
% 
% dg_makeLocalAvgRefs.m Calls dg_localAvgRef for as many local avg refs you
% wish to define.
% 
% dg_subtractLocalAvgRefs.m This does the actual referencing, i.e. it
% subtracts the previously defined local avg ref signals from each of the
% CSC channels that went the avg, and produces a new local avg reffed file
% for each CSC channel.

%$Rev: 239 $
%$Date: 2016-03-10 19:30:16 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $

if nargin < 4
    destdir = sessiondir;
end

% Compute local average refs:
for refidx = 1:length(localavgrefs)
    CSCfilespec = localavgrefs{refidx};
    outfilepath = fullfile(destdir, sprintf( ...
        '%s%d.mat', prefix, refidx ));
    if ~isempty(CSCfilespec)
        dg_localAvgRef(sessiondir, CSCfilespec, outfilepath);
    end
end


