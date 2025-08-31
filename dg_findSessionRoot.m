function sessionroot = dg_findSessionRoot(mydirectory)
%sessionroot = dg_findSessionRoot(mydirectory)
% Finds the root session directory for the session or fragment at
% <mydirectory>.

%$Rev: 146 $
%$Date: 2012-02-06 14:46:53 -0500 (Mon, 06 Feb 2012) $
%$Author: dgibson $

sessionroot = [];
while isempty(sessionroot)
    fraginfofilepath = fullfile(mydirectory, 'lfp_FragmentFiles.mat');
    fraginfofilepath2 = fullfile(mydirectory, 'lfp_fragmentfiles.mat');
    [parent,name] = fileparts(mydirectory);
    if exist(fraginfofilepath, 'file') || exist(fraginfofilepath2, 'file')
        mydirectory = parent;
    else
        if isempty(name)
            % Oops, hit the file system root
            return
        else
            sessionroot = mydirectory;
        end
    end
end
