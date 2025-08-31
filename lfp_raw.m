function lfp_raw(filenum)
%lfp_raw(filenum)

%  Displays the specified CSC file over the entire session with event
%  markers.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if nargin == 0
    filenum =1;
end
lfp_declareGlobals;
fignum = 1000 + filenum;
mytitle = sprintf('%s whole session', lfp_DataDir);
axinfo.xlabel = [];
axinfo.xlim = [];
axinfo.ylim = [];
lfp_multichannel_plot(fignum, mytitle, filenum, ...
    1:numel(lfp_Samples{filenum}), ...
    1:size(lfp_Events,1), 0, axinfo);