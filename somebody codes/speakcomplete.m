clear all
clc;
Fs = 16000;
%sampling frequency is 16000 Hz
[input, fs] = audioread('test.wav');
resampled = resample(input,Fs,fs);
%resample the audio to 16000 Hz
originalLength = length(input)/fs;
soundLength = length(resampled)/Fs;
n = 0.008*Fs;
overlap = 0.006*Fs;

buff= buffer(input,n);
Overlapped = buffer(resampled, n, overlap);

hanningwindow = hann(n);

[rows,cols]=size(Overlapped);

buffandhann= Overlapped;
for j = 1:cols
    for j = 1:rows
        buffandhann(j,j) = Overlapped(j,j)*hanningwindow(j);
    end
end

TH = abs(fft(buffandhann,32));
cutinhalf = TH(1:16,:);
  
  plot(cutinhalf);
  
  testmean = mean(cutinhalf(:,5:20));


%Find the mean of each column, then add all entries of the truncated
%matrix that are larger than the mean to the "peak matrix".
for i = 1:cols
    meanVal = 0.094;
    for j = 1:16
        if (meanVal < cutinhalf(j,i))
            if (cutinhalf(j,i) > 1)
                peakMatrix(j,i) = 1;
           else
                peakMatrix(j,i) = cutinhalf(j,i);
             end
        else
            peakMatrix(j,i) = 0; 
        end
    end
      
end

csvwrite('outputCSV.csv',peakMatrix);
s=vocoder('outputCSV.csv');
sound(s,16000);

% surf(peakMatrix)
% hold on
imagesc(peakMatrix);
