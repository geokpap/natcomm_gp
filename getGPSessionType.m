function [setupstr, sestype] = getGPSessionType(sessiondir)
%OUTPUT
% setupstr: string suitable for submission to 'lfp_changeSetup'.  Empty
%   if setupstr is unrecognized.  <setupstr> for 'ApAv_em' sessions
%   is 'georgios_ApAv'.
% sestype: verbatim session type string extracted from <sessiondir>.
%NOTES
% 20240905: changed name of return value to <setupstr> and added return
%   value <sestype>.  Should be backward-compatible.
% 20241101: added 'ApAv_rev_only'.
% 20241220: added 'ApAp_det_rev'.

setupstr = '';
sessiontypestrings = {
    'ApAp'
    'ApAv'
    'ApAv_em'
    'ApAvApAp'
    'ApAv_rev_only'
    'ApAv_rev_only_ApAv_2bl'
    'ApAp_det_rev'
    };
for typeidx = 1:length(sessiontypestrings)
    sestype = sessiontypestrings{typeidx};
    pat = ['/' sessiontypestrings{typeidx} '/'];
    if ~isempty(strfind(sessiondir, pat)) 
        if isequal(sestype, 'ApAv_em')
            setupstr = 'georgios_ApAv';
        elseif isequal(sestype, 'ApAv_rev_only_ApAv_2bl')
            setupstr = 'georgios_ApAv_rev_only';
        else
            setupstr = ['georgios_' sessiontypestrings{typeidx}];
        end
        return
    end
end
