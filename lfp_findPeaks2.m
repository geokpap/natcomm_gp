function ts = lfp_findPeaks2(fn, thresh1, thresh2)
%ts = lfp_findPeaks2(fn, thresh1, thresh2)  for lfp_createEvents

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

ts = lfp_index2time(lfp_findPeaks(fn, thresh1, thresh2));
