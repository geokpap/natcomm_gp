function dg_MakeMonkeyMatlab(directory)
%dg_MakeMonkeyMatlab(directory)

% Converts the Neuralynx-format monkey data files in directory that are
% needed by lfp_lib into equivalent .MAT files.

%$Rev: 153 $
%$Date: 2012-07-17 18:40:53 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

DataFiles = dir(directory);
for file = DataFiles'
    if ~strcmp(file.name, '.') && ~strcmp(file.name, '..')
        [pathstr,name,ext] = fileparts(file.name);
        if strcmpi(ext, '.NCS') || strcmpi(ext, '.NEV')
            dg_Nlx2Mat(fullfile(directory, file.name));
        end
    end
end
    