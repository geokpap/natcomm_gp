function result = dg_standardizeStringHist(bins, strings, counts)
%DG_STANDARDIZESTRINGHIST converts a histogram of arbitrary string values
%to another histogram that contains a bin for each of a set of standard
%string values.
% <bins> is a list of all string values to be included; <strings> is a list
% of string values actually present; <counts> is the same size as <strings>
% and is the number of instances of each value in <strings>.  <result> is a
% vector of counts that matches <bins>; i.e. if an element of <bins> is not
% present in <strings> it receives a count of 0, and if an element of
% <strings> is not in <bins> it is ignored.  The elements of <strings> are
% assumed to be unique, as are the elements of <bins>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~isempty(strings) && ~isequal(sort(strings), unique(strings))
    error('dg_standardizeStringHist:dupstring', ...
        '<strings> contains at least one duplicated entry.' );
end
if ~isempty(bins) && ~isequal(sort(bins), unique(bins))
    error('dg_standardizeStringHist:dupstring', ...
        '<strings> contains at least one duplicated entry.' );
end
result = zeros(size(bins));
for binnum = 1:length(bins)
    stringidx = find(ismember(strings, bins(binnum)));
    if ~isempty(stringidx)
        result(binnum) = counts(stringidx);
    end
end