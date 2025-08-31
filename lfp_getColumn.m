function result = lfp_getColumn(cellvec, col)
%result = lfp_getColumn(cellvec, col)
%   Gets the contents of one or more columns of a table stored in the
%   format of lfp_EyeTabulation. <cellvec> is a cell vector with one
%   element for each trial.  The cells are assumed to contain matrices with
%   equal numbers of columns, or empty matrices.  If <col> is a number, it
%   is used to index the column in the contents of <cellvec>.  If <col> is
%   a numeric vector, then a column is returned for each element. If <col>
%   is a letter, it is translated to a column number based on the
%   convention that C = 1, D = 2, etc., and A refers to the sequential
%   trial numbers, and B refers to the unique trial IDs.  These are not
%   explicitly stored by lfp_EyeTabulation, and they depend on
%   lfp_SelectedTrials having the same value it did when lfp_EyeTabulation
%   was run.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

result = [];

if strcmp(class(col), 'char')
    if length(col) > 1
        error('lfp_getColumn:badcol', ...
            '<col> must be a single character or a number' );
    else
        col = upper(col);
        if col > 'B'
            col = col - 'B';
        elseif col >= 'A'
            trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
            if length(trials) ~= length(cellvec)
                error('lfp_getColumn:badSelect', ...
                    'Trial selection state has changed since generating <cellvec>' );
            end
            for trialidx = 1:length(cellvec)
                if ~isempty(cellvec{trialidx})
                    if col == 'A'
                        result(end+1,1) = trials(trialidx);
                    else    % col == 'B'
                        result{end+1,1} = lfp_getTrialID(trials(trialidx));
                    end
                end
            end
            return
        end
    end
end

if fix(col) ~= col
    error('lfp_getColumn:badcol2', ...
        '<col> must be a single integer' );
end

nonempty = [];
for trialidx = 1:length(cellvec)
    if ~isempty(cellvec{trialidx})
        nonempty(end+1) = trialidx;
    end
end

myarray = cell2mat(reshape(cellvec(nonempty),[],1));
result = myarray(:, col);
