FOR PHYSIOLOGICAL DATA ANALYSIS (heart rate variability, reaction times, lick activity):
1. addBaselineToMatrix.m  
2. addLickingToMatrix.m  
3. buildMatrixApAv.m  
4. cleanLickByTrial.m  
5. cleanMatrix.m  
6. Code Documentation.docx  
7. loadDataApAv.m  
8. loadLicking.m  
9. matsToCells.m  
10. plotLickingMeans_fullTrial.m


FOR PHYSIOLOGICAL DATA ANALYSIS (pupil diameter):
1. dg_binsearch.m  
2. dg_reconcileFrames.m  
3. dg_substitute.m  
4. dg_zip.m  
5. getGPSessionType.m  
6. lfp_addT1cols3n4.m  
7. lfp_changeSetup.m  
8. lfp_createCSIndices.m  
9. lfp_createTrialIndex.m  
10. lfp_declareGlobals.m  
11. lfp_enabledTrials.m  
12. lfp_findRecSegments.m  
13. lfp_getEvtIDs_georgios_ApAv_rev_only.m  
14. lfp_getEvtIDs.m  
15. lfp_getEyeCalibFileAbs.m  
16. lfp_getSessionInfo.m  
17. lfp_index2time.m  
18. lfp_initializeGlobals.m  
19. lfp_loadSetup.m  
20. lfp_log.m  
21. lfp_otherPreprocess.m  
22. lfp_preamble.m  
23. lfp_read2.m  
24. lfp_readEvents.m  
25. lfp_readEyeCalib.m  
26. lfp_reconcileFrames.m  
27. lfp_selectByRule.m  
28. lfp_setup_georgios_ApAp.m  
29. lfp_setup_georgios_ApAv_rev_only.m  
30. lfp_setup_georgios_ApAv.m  
31. lfp_setup_georgios_ApAvApAp.m  
32. lfp_time2index.m  
33. lfp_trialHasEvent.m  
34. trialByTrialPupilData3.m  


FOR PHYSIOLOGICAL DATA PLOTTING:
1. physiological_data_sim.m  
(loads/plots sample_A.mat and sample_B.mat datasets)

FOR NEURAL DATA:
1. ApAv_GroupSessions_plotWithEventFiles.m  
2. combine_variables.m  
3. ems_segments.m  
4. ft_getopt.m  
5. neural_analysis.m  
6. read_neuralynx_nev.m  
7. readNexFile.m  
8. summaryPlotForApAv_EM.m  


FOR NEURAL DATA PLOTTING:
1. data_analysis.m (loads/plots sample_data.mat dataset)


FOR COMPUTATIONAL MODEL:

1. analysis.py 
2. config.py
3. config_global.py 
4. environment.py
5. logger.py
6. models.py 
7. train.py 
8. train_utils.py 
9. analysis.cpython-39.pyc 
10. config.cpython-39.pyc 
11. config_global.cpython-39.pyc  
12. environment.cpython-39.pyc  
13. logger.cpython-39.pyc 
14. models.cpython-39.pyc  
15. train.cpython-39.pyc  
16. train_utils.cpython-39.pyc  

Neural population-level analysis — ENTRY POINTS (run these):
1. ApAv_GroupSessions_plotWithEventFiles.m
2. summaryPlotForApAv_EM.m
3. neural_analysis.m
4. combine_variables.m
5. ems_segments.m
6. ft_getopt.m
7. read_neuralynx_nev.m
8. readNexFile.m

Neural population-level analysis — FULL helper list (.m)
Because there are hundreds of helper utilities (dg_*.m, lfp_*.m, chronux_*.m),
this TXT ships with a reproducible MATLAB snippet that *appends* the exact list
from your local population-analysis folder (e.g., 'testfart 3') directly into
THIS FILE. Run the following once to finalize the inventory:

%% Append full .m inventory for the population-analysis folder into this TXT.
% Adjust popDir to the folder that contains your dg_*.m / lfp_*.m utilities.
popDir = fullfile(pwd, 'testfart 3');   % <-- change if your folder is named differently
txtFile = fullfile(pwd, 'ALL_FILES_inventory.txt');

d = dir(fullfile(popDir, '*.m'));
names = string({d.name});
names = names(~endsWith(names,"~"));    % drop backup files like *.m~

fid = fopen(txtFile, 'a');              % append to this same TXT
fprintf(fid, '\n-- BEGIN AUTO-GENERATED POPULATION .M INVENTORY (%s) --\n', popDir);
for k = 1:numel(names)
    fprintf(fid, '    %4d. %s\n', k, names(k));
end
fprintf(fid, '-- END AUTO-GENERATED POPULATION .M INVENTORY --\n');
fclose(fid);

fprintf('Appended %d files to %s\n', numel(names), txtFile);

Notes:
- Keep NAS/editor cruft out of Git: @eaDir/, *.m~, *.asv, .DS_Store, Thumbs.db
- Re-run the snippet whenever you add/remove helpers; commit the updated TXT.



The MATLAB scripts are compatible with MATLAB R2018a and newer versions, for both MacOS (Ventura and newer versions) and Windows Operating Systems (Windows 7, 10). Add the MATLAB scripts to your path, load the required data (as specified in the parentheses), and run the scripts. The approximate setup time is 30 minutes.
