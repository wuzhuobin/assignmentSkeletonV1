%codeimplementedfortheformantfreqeuncymethod

clear all; %clears previous work/stored data from workdpace
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

%% %%CIS PROCESSING 

TH = abs(fft(b_h,32));
cutinhalf = TH(1:16,:);
%% %% VOCODER CODING AND ELECTRODOGRAM

csvwrite('outputCSV.csv',cutinhalf);
s=vocoder('outputCSV.csv');
sound(s,Fs);

imagesc(cutinhalf);
