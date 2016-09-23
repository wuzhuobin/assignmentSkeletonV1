function cochlearProc(inName, procType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstrates basic operation of cochlear implant.
% 
% A sound file is read, noise added and then a windowed FFT performed. 
%
% Minimum value filtering and N of M filtering is applied to demonstrate
% the change to the signal. 
%
% note that this demo uses constant chanel width filters. Some (all?) 
% cochlear implants useapproximations of constant Q filters. (Q=f/df)
%
% Version 1.0 - 10-Aug-16
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:   Greg Watkins
% Date:     Jul 2016
% Student:  z5022099
%
% Copyright 2014-6 Greg Watkins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        error('syntax: cochlearProc(inName, procType)')
    end
    
    cochlear = classCochlear();
    
    fprintf('Processing %s with strategy %s\n', inName, cochlear.procName(procType));
    
    % read in speech file
    speech = cochlear.getWav(inName); 
    
    % create a frequency/time matrix
    data = cochlear.getFTM(speech, procType);
    
    % Apply cochlear speech processing
    ftmOut = cochlear.process(data, procType); 
    
    % limit dynamic range etc
    ftmDR  = cochlear.applyDR(ftmOut); 
    
    % plot speech signal
    signalPlot = cochlear.plotSignal(speech); 

    % plot electrodogram
    electrodogram = cochlear.plotElectrodogram(ftmDR, procType); 
    
    % save signal plot to disk
    cochlear.writeJpg( signalPlot, strcat(cochlear.procName(procType), '_Speech'));
    
    % save electrodogram to disk
    cochlear.writeJpg( electrodogram, cochlear.procName(procType) );
 
    % save cochlear stimulation to disk
    cochlear.writeCsv( ftmDR, procType ); 

    % simulate the sound heard after speech processing
    % and save to disk.
    voc    = cochlear.vocoder(ftmDR);
    cochlear.writeWav( voc, procType );
    
end
