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
            if abs(Fs - obj.fTolerance) > 0.1
                y = resample(y, obj.fSample, Fs);
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
            % FTM = zeros()
            if (type == obj.procF0f1f2)
                %
                % extratc formants and build FTM
                % 
            
                % insert your code here
                % How many epochs
                N = floor(length(wav) / obj.fSample / 0.002);
                pointsPerEpoch = length(wav) / N;
                FTM = zeros(obj.numElectrodes, N);
                % Ditch any decimals if length(wav)/N is not an integer!
                pointsPerEpoch = floor(pointsPerEpoch);
                for i = 1:N
                    % get an epoch 
                    epoch(i, :) = wav(((i-1)*pointsPerEpoch) + 1:i*pointsPerEpoch);
                end
                % find the formants
                % the following is a 'rule of thumb' of formant estimation
                numberOfCoefficients = 2 + obj.fSample/1000;
                % parameters initialization
                ffreq = zeros(N, obj.numFormants);
                amplitude = zeros(N);
                for i = 1:N
                    % skip all zeros epoch  
                    if (sum(epoch(i, :)) == 0 )
                        continue;
                    end
                    % determine a polynomial finding the spectral peaks
                    % of the epoch using lpc
                    spectralPeaks = lpc(epoch(i,:)', numberOfCoefficients);
                    r = roots(spectralPeaks); % find roots of this polynomial
                    r = r(imag(r) > 0.01); % only look for + roots up to fs/2
                    % when the r is empty, skip
                    if length(r) < 1
                        continue;
                    end
                    % covert the complex roots to Hz and sort from low to
                    % high
                    formants = sort(atan2(imag(r), real(r))*obj.fSample/(2*pi));
                    % print first three formants (the reset are still stored
                    % in ffreq) 
                    amplitude(i) = max(epoch(i,:));
                    ffreq(i, :) = formants(1: obj.numFormants);
                end
                
                electrodeRange = linspace(obj.fSample/2, 0 ,obj.numElectrodes + 1);
                for i = 1:N
                    for j = 1:obj.numFormants
                        for range = 1:obj.numElectrodes
                            if electrodeRange(range) <= ffreq(i, j) && ffreq(i, j) <= electrodeRange(range + 1)
                                FTM(i, range) = amplitude(i);
                                break;
                            end
                        end
                    end
                end

                    


                
                
            elseif (type == obj.procSpeak) || (type == obj.procCIS)
                %
                % use FFT to implement n x BPF
                %
            
                % insert your code here
                segments = obj.tSample / (1 - obj.frameOverlap) * obj.fSample;
                overlappedSignal = buffer(wav, segments, segments * obj.frameOverlap);
                N = length(overlappedSignal);
                % Hann window filter
%               bsxfun
                hannFilter = repmat(hann(segments), [1, N]);
                % FFT
                frequency = fft(overlappedSignal .* hannFilter, obj.numElectrodes * 2);
                FTM = frequency(1:obj.numElectrodes, :);
            else 
                error('Unknown type (%d)', type)
            end
            
            result = []; % change to return your FTM
            result = FTM;
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