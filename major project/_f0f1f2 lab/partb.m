% Part B&C
close all;clear all;clc;

buffer=sprintf('F1_gregg_hello.wav');
[x,fs]=audioread(buffer);

[n,nChan]=size(x);
if nChan > 1
    error('The type of the wave file must be mono (not stereo)');
end

% epochs
N=100;
points_per_epoch=length(x)/N;
points_per_epoch=floor(points_per_epoch);
for i=1:N
    epoch(i,:)=x(((i-1)*points_per_epoch)+1:i*points_per_epoch);
   
end

% Part D&E
num_of_coefficients = 2+fs/1000;
for i=1:N
    spectral_peaks=lpc(epoch(i,:)',num_of_coefficients);
    r=roots(spectral_peaks);
    r=r(imag(r)>0.01);
    formats=sort(atan2(imag(r),real(r))*fs/(2*pi));
    for j=1:5
        ffreq(i,j)=formats(j);
        fprintf('Epoch(%d), Format(%d)=%.1f Hz\n',i,j,formats(j));
    end
    amplitude(i)=max(epoch(i,:));
end

% Part F
wave1=[];
wave2=[];
wave3=[];
wave4=[];
wave5=[];
for i=1:N
    time=(0:length(epoch(i,:))-1)/fs;
    sound1(i,:)=amplitude(i)*sin(2*pi*ffreq(i,1)*time);
    sound2(i,:)=amplitude(i)*sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time);
    sound3(i,:)=amplitude(i)*sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time);
    sound4(i,:)=amplitude(i)*sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time)+sin(2*pi*ffreq(i,4)*time);
    sound5(i,:)=amplitude(i)*sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time)+sin(2*pi*ffreq(i,4)*time)+sin(2*pi*ffreq(i,5)*time);
    wave1=[wave1 sound1(i,:)];
    wave2=[wave2 sound2(i,:)];
    wave3=[wave3 sound3(i,:)];
    wave4=[wave4 sound4(i,:)];
    wave5=[wave5 sound5(i,:)];
end

% audiowrite('f1.wav',wave1',fs);
% audiowrite('f1f2.wav',wave2',fs);
% audiowrite('f1f2f3.wav',wave3',fs);
% audiowrite('f1f2f3f4.wav',wave4',fs);
% audiowrite('f1f2f3f4f5.wav',wave5',fs);

% Part G
Etordes=4;
stimulus = zeros(N,Etordes);

%% f1
track=1;
for i=1:N
    if ffreq(i,track)>=50 && ffreq(i,track)<=750
        stimulus(i,1)=amplitude(10);
    elseif ffreq(i,track)>=500 && ffreq(i,track)<=1700
        stimulus(i,2)=amplitude(10);
        
    elseif ffreq(i,track)>=1500 && ffreq(i,track)<=4000
        stimulus(i,3)=amplitude(10);
    elseif ffreq(i,track)>=3500 && ffreq(i,track)<=11000 
        stimulus(i,4)=amplitude(10);
    
    end
    
end
figure(1);
imagesc([1,N],[1:Etordes],(stimulus'));
colorbar


%% f2
stimulus = zeros(N,Etordes);
track=2;
for i=1:N
    if ffreq(i,track)>=50 && ffreq(i,track)<=750
        stimulus(i,1)=amplitude(10);
    elseif ffreq(i,track)>=500 && ffreq(i,track)<=1700
        stimulus(i,2)=amplitude(10);
        
    elseif ffreq(i,track)>=1500 && ffreq(i,track)<=4000
        stimulus(i,3)=amplitude(10);
    elseif ffreq(i,track)>=3500 && ffreq(i,track)<=11000 
        stimulus(i,4)=amplitude(10);
    
    end
    
end
figure(2);
imagesc([1,N],[1:Etordes],(stimulus'));
colorbar