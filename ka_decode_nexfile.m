function [nexFile,neuronCount,eventCount,intervalCount,waveCount,popCount,contCount,markerCount]=ka_decode_nexfile(openfilename0)
%
%NOTES
% Received from Ken-Ichi Amemori 22-Oct-2015, nothing changed but the name
% and header comments (it was originally just "decode_nexfile").

%$Rev: 231 $
%$Date: 2015-10-22 18:55:01 -0400 (Thu, 22 Oct 2015) $
%$Author: dgibson $

% calc max trial num
fileName=openfilename0
nexFile = [];

fid = fopen(fileName, 'r');
if(fid == -1)
    error 'Unable to open file'
    save ('temp.mat','stats')
    clear variables;
    load temp.mat
    delete temp.mat
%   clear global -regexp stats;
    return
end

magic = fread(fid, 1, 'int32');
if magic ~= 827868494,
    error 'The file is not a valid .nex file'
end;
nexFile.version = fread(fid, 1, 'int32');
nexFile.comment = deblank(char(fread(fid, 256, 'char')'));
nexFile.freq = fread(fid, 1, 'double');
nexFile.tbeg = fread(fid, 1, 'int32')./nexFile.freq;
nexFile.tend = fread(fid, 1, 'int32')./nexFile.freq;
nvar = fread(fid, 1, 'int32');

% skip location of next header and padding
fseek(fid, 260, 'cof');

neuronCount = 0;
eventCount = 0;
intervalCount = 0;
waveCount = 0;
popCount = 0;
contCount = 0;
markerCount = 0;

% real all variables
for i=1:nvar
    type = fread(fid, 1, 'int32');
    varVersion = fread(fid, 1, 'int32');
	name = deblank(char(fread(fid, 64, 'char')'));
    offset = fread(fid, 1, 'int32');
	n = fread(fid, 1, 'int32');
    WireNumber = fread(fid, 1, 'int32');
	UnitNumber = fread(fid, 1, 'int32');
	Gain = fread(fid, 1, 'int32');
	Filter = fread(fid, 1, 'int32');
	XPos = fread(fid, 1, 'double');
	YPos = fread(fid, 1, 'double');
	WFrequency = fread(fid, 1, 'double'); % wf sampling fr.
	ADtoMV  = fread(fid, 1, 'double'); % coeff to convert from AD values to Millivolts.
	NPointsWave = fread(fid, 1, 'int32'); % number of points in each wave
	NMarkers = fread(fid, 1, 'int32'); % how many values are associated with each marker
	MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
	MVOfffset = fread(fid, 1, 'double'); % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
    filePosition = ftell(fid);
    
    switch type
        case 0 % neuron
            neuronCount = neuronCount+1;
            nexFile.neurons{neuronCount,1}.name = name;
            fseek(fid, offset, 'bof');
            nexFile.neurons{neuronCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            fseek(fid, filePosition, 'bof');
            
        case 1 % event
            eventCount = eventCount+1;
            nexFile.events{eventCount,1}.name = name;
            fseek(fid, offset, 'bof');
            nexFile.events{eventCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            fseek(fid, filePosition, 'bof');
        
        case 2 % interval
            intervalCount = intervalCount+1;
            nexFile.intervals{intervalCount,1}.name = name;
            fseek(fid, offset, 'bof');
            nexFile.intervals{intervalCount,1}.intStarts = fread(fid, [n 1], 'int32')./nexFile.freq;
            nexFile.intervals{intervalCount,1}.intEnds = fread(fid, [n 1], 'int32')./nexFile.freq;
            fseek(fid, filePosition, 'bof');  
        
        case 3 % waveform
            waveCount = waveCount+1;
            nexFile.waves{waveCount,1}.name = name;
            nexFile.waves{waveCount,1}.NPointsWave = NPointsWave;
            nexFile.waves{waveCount,1}.WFrequency = WFrequency;
            
            fseek(fid, offset, 'bof');
            nexFile.waves{waveCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            nexFile.waves{waveCount,1}.waveforms = fread(fid, [NPointsWave n], 'int16').*ADtoMV + MVOfffset;
            fseek(fid, filePosition, 'bof'); 
            
        case 4 % population vector
            popCount = popCount+1;
            nexFile.popvectors{popCount,1}.name = name;
            fseek(fid, offset, 'bof');
            nexFile.popvectors{popCount,1}.weights = fread(fid, [n 1], 'double');
            fseek(fid, filePosition, 'bof');
            
        case 5 % continuous variable
            contCount = contCount+1;
            nexFile.contvars{contCount,1}.name = name;
            nexFile.contvars{contCount,1}.ADFrequency = WFrequency;
            fseek(fid, offset, 'bof');
            nexFile.contvars{contCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            nexFile.contvars{contCount,1}.fragmentStarts = fread(fid, [n 1], 'int32') + 1;
%            nexFile.contvars{contCount,1}.data = fread(fid, [NPointsWave 1], 'int16').*ADtoMV + MVOfffset;
            fseek(fid, filePosition, 'bof'); 
         
        case 6 % marker
            markerCount = markerCount+1;
            nexFile.markers{markerCount,1}.name = name;
            fseek(fid, offset, 'bof');
            nexFile.markers{markerCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            for itmp=1:NMarkers
                nexFile.markers{markerCount,1}.values{itmp,1}.name = deblank(char(fread(fid, 64, 'char')'));
                for p = 1:n
                    nexFile.markers{markerCount,1}.values{itmp,1}.strings{p, 1} = deblank(char(fread(fid, MarkerLength, 'char')'));
                end
            end;
            fseek(fid, filePosition, 'bof');
        
        otherwise
            disp (['unknown variable type ' num2str(type)]);
    end
	dummy = fread(fid, 60, 'char');
end
fclose(fid);

