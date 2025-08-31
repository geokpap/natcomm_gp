function dg_CreateMVSpreadsheet(directory, files, segspecs, segspecstrings, ...
    condition, selectionspec)

%dg_CreateMVSpreadsheet(directory, files, segspecs, segspecstrings, condition)
% Creates an Excel-readable spreadsheet containing statistical measures for
% each segment in segspecs for each unit in each file in files located in
% directory directory.  files is a cell array of strings.  segspecs is a
% vector of segment specification structures of the type returned by
% dg_ParseSegmentSpec.  segspecstrings is a cell string array of the
% segment specifications.
% condition is the thing returned by dg_ParseSelectionSpec(selectionspec).
% The 'Measure' column in the spreadsheet indicates which statistic is
% reported on each row; 'normvar' is an abbreviation for 'normalized
% variance', which is here defined as variance/(mean^2).

%$Rev: 153 $
%$Date: 2012-07-17 18:40:53 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

[OutputFileName, OutputPathName] = uiputfile(...
    fullfile(directory, 'MeanVar.XLS') );
outfile = fopen(fullfile(OutputPathName, OutputFileName), 'w');
fprintf(outfile, 'Selection Spec = "%s"\n', selectionspec);
fprintf(outfile, 'Filename\tSession\tUnit ID\tMeasure');
for specstring = segspecstrings'
    fprintf(outfile, '\t%s', char(specstring));
end
fprintf(outfile, '\n');
measurelabels = [ {'mean'} {'variance'} {'normvar'}];
for file = files'
    fullfilename = fullfile(directory, char(file));
    [fileheader, trialdata] = dg_ReadMouseFormat(fullfilename);
    [pathstr, session, electrode] = fileparts(fullfilename);
    session = lower(session);
    electrode = lower(electrode);
    electrode = electrode(3:4);
    if electrode(1) == 't' || electrode(1) == '0'
        electrode = electrode(2);
    end
    for clusternum = (1:fileheader.CSize)
        rates = dg_RatePerTrial(segspecs, trialdata, clusternum, condition);
        MeanAndVariance = dg_MeanAndVariance(rates);
        MeanAndVariance = [ MeanAndVariance
            MeanAndVariance(2,:)./(MeanAndVariance(1,:).^2) ];
        for r = (1:3)
            fprintf(outfile, '%s\t%s\tT%sU%d\t%s', fullfilename, ...
                session, electrode, clusternum, char(measurelabels(r)));
            for k = (1:length(MeanAndVariance(r,:)))
                fprintf(outfile, '\t%.2E', MeanAndVariance(r,k));
            end
            fprintf(outfile, '\n');
        end
    end
end
fclose(outfile);
