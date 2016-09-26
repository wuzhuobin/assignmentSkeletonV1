%formant frequency
clear all;
clc;

%sread the input sound and store data in in_s, and the sampling freqeuncy
%in fs
[in_s, fs] = audioread('test.wav');
%resample the sound into a sampling freqeuncy of 16000. So this program
%will convert any sampling frequency of any input audio to 16000Hz.
Fs = 16000;
a = resample(in_s,Fs,fs);

originalLength = length(in_s)/fs;
soundLength = length(a)/Fs;
n = 0.008*Fs;
overlap = 0.006*Fs;

buff= buffer(in_s,n);
Overlapped = buffer(a, n, overlap);

hanningwindow = hann(n);

[rows,cols]=size(Overlapped);

buffandhann= Overlapped;
for i = 1:cols
    for j = 1:rows
        buffandhann(j,i) = Overlapped(j,i)*hanningwindow(j);
    end
end

TH = abs(fft(buffandhann,32));
cutinhalf = TH(1:16,:);

csvwrite('outputCSV.csv',cutinhalf);
s=vocoder('outputCSV.csv');
sound(s,Fs);

% mesh(cutinhalf); % interpolated
% axis tight; hold on
% imagesc(cutinhalf)
