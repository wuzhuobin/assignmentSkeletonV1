classdef classCochlearSupport < handle 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This class provides data definitions and methods to support the 
% implementation of a cochlear implant sound porcessor simulation.
%
% Methods:
%   procName
%   vocoder
%   writeJpg
%   writeWav
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

    properties
    end
    
    properties (Constant)
    %
    % contstant data required for cochlear implant
    %
    % data can be access as obj.<name> eg obj.fs
    %
        procF0f1f2   = 1; 
        procSpeak    = 2;
        procCIS      = 3;
        procNames    = {'F0f1f2';
                        'Speak';
                        'CIS';};
        
        fSample      = 16000; % required sampling frequency
        fTolerance   = 0.1;   % allowed variation from fSample
        
        tSample      = 0.002; % sample time for cochlear implant processing
        frameOverlap = 3/4;   % required overlap when using Hann window
        
        numElectrodes= 16;    % number of electrodes in cochlear implant
        numFormants  = 2;     % formants used by the F0F1F2 strategy
                              
        dynamicRange = 10;    % dB
        maxOutput    = 1024;
        threshold    = 0.1;
        
        dpi          = '-r300'; % dots per inch resolution for saved image files
        imageType    = '-djpeg';% type of graphics file for saved images.
        jpgType      = '.jpg';  % 
        wavType      = '.wav';  % 
        csvType      = '.csv';  %
        
        fontSize     = 14;      % good size fore labelling electrodograms

    end
    
    properties % (Access = private)

    end

    methods 

        %--------------------------------------------------------
        % initialise class
        %--------------------------------------------------------
        function obj = classCochlearSupport()
            % all initializations, calls to base class, etc. here,
        end
        

        %--------------------------------------------------------
        % methods
        %--------------------------------------------------------

        function result = procName(obj, procType)
        %
        % return string containing the name of the process type 
        % corresponding to <procType> 
        %
            if (procType < obj.procF0f1f2) || (procType > obj.procCIS)
                error('Unknown type (%d)', procType)
            end
            result = obj.procNames{procType};
        end
        
        function result = vocoder(obj, ftm)
        %
        % create a vocoded output from the <ftm>
        %
        % the sample time per column is tSample. The frame overlap % is
        % frameOverlap
        %
        % return the wav vector
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);
            result = vocoder(ftm, obj.tSample, obj.fSample);
        end
        
        function result = writeJpg(obj, h, name )
        %
        % prints the figure with handle <h> to a jpg file 
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);
            name = strcat(pwd, '\', name, obj.jpgType);
            print(h, name, obj.imageType, obj.dpi);
            result = [];
        end
        
        function result = writeWav(obj, data, procType)
        %
        % writes the data to a wav file
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);
            name = strcat(pwd, '\', obj.procName(procType), obj.wavType);
            if isempty(data)
                warning('wav data was empty');
            else
                audiowrite(name, data, obj.fSample);
            end
            result = [];
        end

        function result = writeCsv(obj, data, procType )
        %
        % writes the data to a csv file
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);
            name = strcat(pwd, '\', obj.procName(procType), obj.csvType);
            csvwrite(name, data);
            result = [];
        end

    end
end

function [result] = vocoder(freqSamples, sampleLength, fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vocoder recreates signal as it might be heard by a CI recipient.
%
% INPUT:
%   freqSamples:    The input paramter is either a Frequency Time Matrix 
%                   (FTM - see below) or the name of a CSV file containing
%                   an FTM. The filename may be a full path including file 
%                   name or just the name of a file that existsin in the
%                   path.
%
%                   The FTM contains one row per electrode. One column is 
%                   present for each time sample.
%
%                   The duration of each time slice (column) and the sample
%                   frequency must be defined for the vocoder. These 
%                   settings are configured as below:-
%                       sampleLength: Duration of each time slice in seconds
%                       fs          : Rate in Hz at which the sound was
%                                      originally sampled.
%                   The channels are assumed to have the same bandwidth.
%                   Notes:-
%                    1. The first channel coresponds to 0 Hz! 
%                    2. In real cochlear implants the channels do not
%                       necessarilly have the same bandwidth.
%
% OUTPUT:
%   result:         vocoded approximation of the sound represneted by the
%                   FTM. The sample rate of the data is fs. This vector can
%                   be played with the command "sound(result,16e3)"
%
% DESCRIPTION
% The vocoder multiplies each column of data by a column of carriers and
% then sums the column to have a waveform that apporximates the sound in
% that time slice.
%
% This operation is repeated for all columns and the data concatenated.
%
% see http://sethares.engr.wisc.edu/vocoders/channelvocoder.html
%
% CARRIER
% Three types of carrier are supported:-
%       carrierSine        : pure sine wave
%       carrierNoiseLinear : linearly spaced set of sine waves effectively
%                            providing full band noise
%       carrierNoiseRandom = randomly distributed set of sine waves 
%                            effectively providing full band noise
% Which carrier is best?.... seems to be subjective.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:   Greg Watkins
% Date:     Aug 2015
% Student:  z5022099
%
%   ed 1    : 30/8/15
%
% Copyright 2014-5 Greg Watkins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if ischar(freqSamples)
%         freqSamples = csvread(freqSamples);
%     end
%     
%     sampleLength = 0.002;
%     fs           = 16e3;
% 
    bw       = fs/2;
    channels = size(freqSamples,1);
    
    carrierSine         = true;
    carrierNoiseLinear  = false;
    carrierNoiseRandom  = false;

    %
    % Filter/carrier information
    % 
    channelWidth = bw/channels;
    centreF = 0:channelWidth:bw-1;
    centreF = centreF';

    tc = 0:1/fs:sampleLength-1/fs;
    
    %
    % Create acrriers
    %
    if carrierSine == true
        carrier = sin(2*pi*centreF*tc);
    elseif carrierNoiseLinear == true
        carrier = zeros(size(centreF,1), size(tc,2));
        n = 1;
        for f=0:channelWidth:fs/2-1
            fc = linspace(f-channelWidth*.4, f+channelWidth*.4,10)';
            tmp = sum(sin(2*pi*fc*tc));
            tmp = tmp/max(tmp);
            carrier(n,:) = tmp;
            n = n + 1;
        end
    elseif carrierNoiseRandom == true
        carrier = zeros(size(centreF,1), size(tc,2));
        n = 1;
        for f=0:channelWidth:fs/2-1
            fc = rand(10,1).*channelWidth-channelWidth/2 + f;
            tmp = sum(sin(2*pi*fc*tc));
            tmp = tmp/max(tmp);
            carrier(n,:) = tmp;
            n = n + 1;
        end
    else
        error('Unknown carrier type');
    end
    
    %
    % synthesise signal
    %
    voc = [];
    for n=1:size(freqSamples,2)
        sample = zeros(size(carrier));
        for m=1:size(freqSamples,1)
            sample(m,:) = carrier(m,:)*freqSamples(m,n);
        end
        voc = [voc sum(sample)];
    end
    
    result = voc;
    result = result';
    result = result/max(result);  % normalise
    result = result*0.999; % as wav file max value must be < 1
end