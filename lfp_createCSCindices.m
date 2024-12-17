function lfp_createCSCindices

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

[startSampleIndex, endSampleIndex] = ...
    lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];

[startrecsample, endrecsample, recOnSamples, recOffSamples] = lfp_findRecSegments;
if ( isempty(startrecsample) ...
        || isempty(endrecsample) ...
        )
    error('lfp_createCSCindices:recsamples');
end
lfp_TrialRec = [ startrecsample endrecsample ];
lfp_RecSegments = [ recOnSamples recOffSamples ];

