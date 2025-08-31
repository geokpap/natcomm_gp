function [startidx, endidx] = dg_findruns(v)
%INPUTS
% v: logical row vector.
%OUTPUTS
% startidx: the index into <v> of the first element in each continuous run
%   of consecutive <true> values.
% endidx: the index into <v> of the last element in each continuous run
%   of consecutive <true> values.  Same length as <startidx>.
%NOTES
% <startidx> and <endidx> can contain any values in the range 1:length(v).
% A "run" of just a single consecutive <true> value is still considered to
% be a "run", in which case <startidx> and <endidx> will be the same for
% that run.

%$Rev: 305 $
%$Date: 2023-09-08 14:19:01 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

startidx = find([ v(1) v(2:end) & ~v(1:end-1) ]);
endidx = find(~v(2:end) & v(1:end-1));
% endidx(1) indexes the last point in the first run, and is necessarily >=
% startidx(1).
if length(startidx) > length(endidx)
    % the last run goes to the end of the <v>.
    endidx(end+1) = length(v);
end

