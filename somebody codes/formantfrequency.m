function [] = formantfrequency(audio)
%codeimplementedforFORMANTFREQUENCYmethod

%clear all; %clears previous work/stored data from workdpace
clc;%clears command window

[X, fs] = audioread('test.wav');

%% PREPROCESSING
Fs = 16000; % resamplingfreqquency

resampled = resample(X,Fs,fs);
%resample the audio to 16000 Hz
o_length = length(X)/fs; % this is just measuring the length of the original file in seconds
RS_length = length(resampled)/Fs;% this measures the length of the resampled sound in seconds
segments = 0.008*Fs;% 0.008 represents 8miliseconds so we are taking 8 milimeter segments of the resampled audio
%3ms before and after audio were requested in this project as an overlap to the 2ms
%which is why the sm2 segement is taken as 8ms in the previous line (3overlap +2msoriginasegment+ 3overlap) 
Requested_overlap = 0.006*Fs; 

%y = buffer(x,n,p) overlaps or underlaps successive frames in the output matrix by p samples:
Overlapped_singal = buffer(resampled, segments, Requested_overlap);

%The Hann function is used as a window function in digital signal 
%processing to select a subset of a series of samples to perform 
%a Fourier transform or other calculations, w = hann(L) returns an L-point symmetric Hann
%window in the column vector w. L must be a positive integer.
H_W = hann(segments);

%create a matrix as the sze of the matrix you intend to process,
%(preserving dimensions)

[rows,cols]=size(Overlapped_singal);
%Apply both buffer and Hann to the signal with overlapping data. 
b_h= Overlapped_singal;
for i = 1:cols
    for j = 1:rows
        b_h(j,i) = Overlapped_singal(j,i)*H_W(j);
    end
end
%% FORMANT FREQUENCY DETECTION
%createnewmatrixwith no columns conatining all zeroes
count = 1;
cols = size(b_h,2);
for c = 1:cols
   if sum(b_h(:,c)) ~= 0
      no_zero(:,count) = b_h(:,c);
      count = count + 1;
   end
end
%calculate the formants
[R,C] = size(no_zero);
number_of_coefficients = 2+fs/1000;
for i=1:C
    % determine a polynomial finding the spectral peaks
    % of the epoch using lpc
    spectral_peaks = lpc(no_zero(:,i)',number_of_coefficients);
    r=roots(spectral_peaks); % find roots of this polynomial
    r=r(imag(r)>0.01); % only look for + roots up to fs/2
    % convert the complex roots to Hz and sort from low to high.
    formants=sort(atan2(imag(r),real(r))*fs/(2*pi));
    % print first five formants (the rest are still stored in ffreq)
    for j=1:5
        ffreq(i,j) = formants(j); % save for later
        
        %fprintf('Epoch(%d), Formant(%d) =%.1f Hz\n',i,j,formants(j));
    end
    amp(i) = max(no_zero(:,i));
end

%% ASSIGNING ELECTRODES 
% only using 8 electrodes because max frequency is 8000 Hz, according to
% the Nyquist sampling rate concept. Nyq maxFrequency = Fs/2
electrod_stimuli = zeros(8,C);
frequencies = [0:1000:8000];
for i = 1:C
    for j = 1:5
       for f = 1:8
          if (ffreq(i,j) > frequencies(f)) && (ffreq(i,j) <= frequencies(f+1))
              electrod_stimuli(f,i) = amp(i);
          end
       end
    end
end
%% %% WRITE THE FILE TO THE VOCODER FORMAT AND PLOT THE ELECTRODOGRAM

csvwrite('outputCSV.csv',electrod_stimuli);
s=vocoder('outputCSV.csv');
sound(s,16000) %playsthesound

imagesc(electrod_stimuli);
