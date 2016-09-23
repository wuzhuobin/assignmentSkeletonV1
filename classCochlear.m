classdef classCochlear < handle & classCochlearSupport 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This class implements methods for a cochlear implant sound processor 
% simulation.
%
% Methods:
%   getWav
%   getFTM
%   process
%   applyDR
%   plotSignal
%   plotElectrodogram
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

    methods 

        %--------------------------------------------------------
        % initialise class
        %--------------------------------------------------------
        function obj = classCochlear()
            % all initializations, calls to base class, etc. here,
        end
        
        
        %--------------------------------------------------------
        % functional methods
        %--------------------------------------------------------

        function result = getWav(obj, name)
        %
        % read the wav file <name> 
        %
        % resample if required so that the sampling frequency is 
        % fs +/- fTolerance
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            % insert your code here
            result = [];
            [y,Fs] = audioread(name);
            % resample if required so that the sampling frequency is 
            % fs +/- fTolerance
            if abs(Fs - fTolerance) > 0.1
                y = resample(y, fSample, Fs);
            end
            result = y;
        end
        
        function result = getFTM(obj, wav, type)
        %
        % return a Frequency Time matrix where each row maps to an
        % electrode and each column is a sample taken at tSample.
        %
        % The apical electrode is in row 1.
        %
        % A Hann window must be applied before an FFT to avoid spurious 
        % high frequency artefacts.
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            % insert your code here
            % Hann window filter
            hannFilter = hann(fSample);
            
            if (type == obj.procF0f1f2)
                %
                % extratc formants and build FTM
                % 
            
                % insert your code here
                

                % FFT
                frequency = fft(wav .* hannFilter);
                
                
            elseif (type == obj.procSpeak) || (type == obj.procCIS)
                %
                % use FFT to implement n x BPF
                %
            
                % insert your code here
            
            else 
                error('Unknown type (%d)', type)
            end
            
            result = []; % change to return your FTM
        end
        
        function result = process(obj, data, type)
        %
        % Apply cochlear sound processing as defined by <type> to <data> 
        % and return the modified data.
        %
        % <data> will be a wav vector for formant based processing; 
        % otherwise an ftm.
        %
        % The different processing types are defined by procXxxxx in the
        % constant sections of the class.
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            if (type == obj.procF0f1f2)
                %
                % Implement F0f1f2 processing
                % 
            
                % insert your code here
            
            elseif (type == obj.procSpeak)
                %
                % Implement SPEAK processing.
                %

                % insert your code here
                
            elseif (type == obj.procCIS)
                %
                % Implement CIS processing.
                %

            
                % insert your code here
            
            else 
                error('Unknown type (%d)', type)
            end

            result = []; % change to return your FTM
        end
        
        function result = applyDR(obj, ftm)
        %
        % apply dynamic rage limiting to dynamicRange. Optionally other
        % scaling can be applied if intelligibility is improved.
        %
        % return the modified ftm.
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            % insert your code here

            result = []; % change to return your FTM
        end
        
        function result = plotSignal(obj, speech)
        %
        % create plot of speech signal. 
        % label axes and apply title.
        %
        % return the figure handle
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            result  = figure;
            
            % insert your code here
            
        end
        function result = plotElectrodogram(obj, ftm, procType)
        %
        % cretae plot of electrodogram with apical electrode at top of
        % plot. <type> defines the type of processing that has been applied.
        % label axes and apply title.
        %
        % return the figure handle
        %
            db=dbstack(); fprintf('    >>%s\n', db(1).name);

            result = figure;
            
            % insert your code here
            
        end
    end
end