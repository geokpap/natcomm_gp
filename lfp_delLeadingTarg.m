function s2 = lfp_delLeadingTarg(s1, targnum)
%   If <targnum> contains more than one element, then any leading targets
%   that are NOT members of targnum are deleted.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

s2 = cell(size(s1));

for trialidx = 1:length(s1)
    if ~isempty(s1{trialidx})
        fixations = s1{trialidx}(:,2);
        if isequal(size(targnum), [1 1])
            nontargnum = find(fixations ~= targnum);
            if isempty(nontargnum)
                nontargnum = size(fixations);
            end
            remove = find(fixations(1:nontargnum(1)) == targnum);
        else
            istargnum = find(ismember(fixations, targnum));
            if isempty(istargnum)
                istargnum = size(fixations);
            end
            remove = find(~ismember(fixations(1:istargnum(1)), targnum));
        end
        s2{trialidx} = s1{trialidx};
        s2{trialidx}(remove,:) = [];
    end
end