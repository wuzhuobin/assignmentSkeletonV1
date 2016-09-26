close all;clear all;clc;
% singal frequency
frequency = 4000;   % Hz

duration = 2;       % second
fs = 44000;            % sampling f Hz
time = [0:(duration*fs)-1];
x = cos(frequency*2*pi*time/fs);
buffer=sprintf('TIME DOMAIN PLOT OF PURE FREQUENCY SIGNAL');
fig1=figure(1);
plot(time/fs,x,'-x');
title(buffer);
axis([0 0.005 -1 1]);
xlabel('Time (s)');
ylabel('Amplitude');
% print(fig1,'my_plot','-dpng');

% audiowrite('my_wave.wav',x',fs);

fft_x=abs(fft(x,length(x)));

frequency_range=fs*[0:length(x)-1]/length(x);
figure(2);
plot(frequency_range,fft_x);
xlabel('Frequency (Hz)');
ylabel('FFT Magnitidue');
title('FREQUENCY DOMAIN PLOT OF PURE FREQUENCY SIGNAL')
axis([0 fs/2 0 max(fft_x)]);



