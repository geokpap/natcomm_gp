function midTS = lfp_markMidMask(mask, varargin)
%midTS = lfp_markMidMask(mask)
% Returns time stamps of the middles of masked episodes, which are defined
% to be consecutive runs of the value <true>.
%INPUTS
% mask:  a logical array
%OUTPUTS
% midTS:  a column vector of timestamps
% OPTIONS:
% 'medplus' - calculates the median duration of the masked episodes and
%   invokes 'mindur' with <mindur> set to the median duration.
% 'mindur', mindur - only include masked episodes whose duration is at
%   least <mindur>.

%$Rev: 236 $
%$Date: 2011-06-11 19:00:04 -0400 (Sat, 11 Jun 2011) $
%$Author: dgibson $

lfp_declareGlobals;

medplusflag = false;
mindur = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'medplus'
            medplusflag = true;
        case 'mindur'
            argnum = argnum + 1;
            mindur = varargin{argnum};
        otherwise
            error('lfp_markMidMask:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

isonset = mask(2:end) & ~mask(1:end-1);
isonset = [ false reshape(isonset, 1, []) ];
isoffset = ~mask(2:end) & mask(1:end-1);
isoffset = [reshape(isoffset, 1, []) false ];

onsamps = find(isonset);
offsamps = find(isoffset);
pairs = dg_zip(onsamps, offsamps);
episodebounds = [ reshape(onsamps(pairs(:,1)), [], 1) ...
    reshape(offsamps(pairs(:,2)), [], 1) ];
midTS = mean(lfp_index2time(episodebounds),2);

if medplusflag
    mindur = median(lfp_index2time(episodebounds(:,2)) - ...
        lfp_index2time(episodebounds(:,1)));
end

if mindur > 0
    longenuff = lfp_index2time(episodebounds(:,2)) - ...
        lfp_index2time(episodebounds(:,1)) >= mindur;
    midTS = midTS(longenuff);
end

