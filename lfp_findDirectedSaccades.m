function saccstarts = lfp_findDirectedSaccades(evtbounds, anglebounds, ...
    saccades)
%saccstarts = lfp_findDirectedSaccades(bounds, saccades)
%   Requires that a value be assigned to the variable 'SaccStart' in your
%   lfp_getEvtIDs_* file.
% INPUTS
%   evtbounds: a pair of indices into lfp_Events
%   anglebounds: lower (first element) and upper (second element) bounds on
%   saccade direction, in degrees relative to X-axis (negative values are
%   automatically converted to corresponding positive values).  If the
%   upper bound is less than the lower bound, then the allowed range of
%   directions includes 0 degrees.
% OUTPUTS
%   saccstarts: a list of saccade start times that took place between the two
%   events in <evtbounds>, and whose directions were in the closed interval
%   specified by <anglebounds>.  Saccades whose direction was NaN are not
%   included.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_Events;
lfp_getEvtIDs;
if exist('SaccStart') ~= 1
    error('lfp_findDirectedSaccades:noSaccStart', ...
        'There is no SaccStart variable defined in your lfp_getEvtIDs_* file.' );
end
if isempty(SaccStart) 
    error('lfp_findDirectedSaccades:noSaccStart2', ...
        'SaccStart must have a non-empty value.' );
end

anglebounds = mod(anglebounds, 360);
if anglebounds(1) > anglebounds(2)
    angleboundses = {[0 anglebounds(2)] [anglebounds(1) 360]};
else
    angleboundses = {anglebounds};
end

timebounds = lfp_Events(evtbounds, 1);

saccadetable = cell2mat(reshape(saccades, [], 1));
saccidx = find(saccadetable(:,4) > timebounds(1) & ...
    saccadetable(:,4) < timebounds(2) );
saccidx2 = [];
angles = mod(saccadetable(saccidx,1), 360);
for k = 1:length(angleboundses)
    saccidx2 = [ saccidx2
        saccidx(find( angles >= angleboundses{k}(1) ...
        & angles <= angleboundses{k}(2) )) ];
end
saccstarts = saccadetable(saccidx2, 4);
if length(angleboundses) > 1
    saccstarts = sort(saccstarts);
end

