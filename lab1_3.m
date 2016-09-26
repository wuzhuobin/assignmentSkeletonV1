clear
% the filename to process
buffer = sprintf('Hello.wav');
% find out sampling rate, etc
[x, fs] = audioread(buffer); %reads the data from the buffer file and return
                            %sampled data x, and a sample rate for that data, fs
% x = waveform vector;
% fs = sampling frequency;
% Make sure that the file is in the correct format before proceeding
[n, nChan] = size(x) % returns n rows and nChan columns if x is a matrix
if nChan > 1
        error ('The type of the wave file must be mono (not stereo)');
end

% plot this sample
duration = length(x);
time = [0:(duration)-1];
%subplot(m,n,p) divides the current figure into an m-by-n grid and 
%creates an axes for a subplot in the position specified by p.
subplot (2,1,1),
plot(time/fs,x);
% time here is the number of data points sampled,
% x here has the unit of number of data points per sec
% hence time/fs has the unit of sec (time)

% plot the second one, frequencies against magnitude
fft_x = abs(fft(x, length(x)));
frequency_range = fs*[0:length(x)-1]/length(x);
subplot(2,1,2),plot(frequency_range, fft_x);
axis([0 fs/2 0 max(fft_x)]);

% the part c of this experiment 

% the signal duration is 2 minutes
% Choose 20 ms to be our epoch
% So we will have 100 epochs

% How many epochs?
N=100;
% How many datapoints in each epoch?
points_per_epoch = length(x)/N;
% Ditch any decimals if length(x)/N is not an integer!
points_per_epoch = floor(points_per_epoch);
for i = 1:N
        % get an epoch
        epoch (i,:) = x(((i-1)*points_per_epoch)+1:i*points_per_epoch);
end

% Part D: Formants of the epochs

% Find the formants
% The following is a 'rule of thumb' of formant estimation
number_of_coefficients = 2+fs/1000;
for i=1:N
    % determine a polynomial finding the spectral peaks of the epoch using
    % lpc
    spectral_peaks = lpc(epoch(i,:)',number_of_coefficients);
    amplitude(i) = max(epoch(i,:)); 
    % find roots of this polynomial
    r = roots(spectral_peaks);
    r = r(imag(r)>0.01); %only look for + roots up to fs/2
                         %imag returns only the imaginary part of the r
    % Convert the complex roots to Hz and sort from low to high.
    % sort (x,dim) returns the sorted elements of A along dimension dim
    % atan(x) returns the inverse tangent of the element of x
    % atan2(x,y) converts variables from cartesian coordinate to polar
    % coordinate, where phi=atan(y,x)
    formants = sort(atan2(imag(r),real(r))*fs/(2*pi));
    % print first five formants )the rest are still stored in ffreq)
    for j=1:5
        ffreq(i,j) = formants (j); %solve for later use
    % ffreq is a matrix, formants for each epoch are stored in this matrix
    %    fprintf('Epoch(%d),Formant(%d)=%.1 Hz\n',i,j,formants(j));
    % the above line actually works, but just for avoding the values from
    % being printed, a % symbol is added in front of the line
        % fprintf function prints multiple numeric values and literal text
        % to the screen
        % Formant (%d)= %.1 Hz is in the formatspec input, prints each
        % value in the vector
        %\n is a control character that starts a new line
    end
end

% part F: Building the sounds
% process the sound based on the formant of each epoch
for i=1:N
        % time for each sample
        time = (0:length(epoch(i,:))-1)/fs; % sampling time
        % Construct a  signal using f1 only to  represent this epoch
        sound1(i,:) = amplitude(i) * sin(2*pi*ffreq(i,1)*time);
        sound2(i,:) = amplitude(i) * (sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time));
        sound3(i,:) = amplitude(i) * (sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time));
        sound4(i,:) = amplitude(i) * (sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time)+sin(2*pi*ffreq(i,4)*time));
        sound5(i,:) = amplitude(i) * (sin(2*pi*ffreq(i,1)*time)+sin(2*pi*ffreq(i,2)*time)+sin(2*pi*ffreq(i,3)*time)+sin(2*pi*ffreq(i,4)*time)+sin(2*pi*ffreq(i,5)*time));
end

% Add all the sounds together in series to build a new sound, 2s long
% Based only on the pure formant frequencies

% Create an empty matrix that we'll add to below
wave1 = 0;
wave2 = 0;
wave3 = 0;
wave4 = 0;
wave5 = 0;
% piece together all the epochs into one sound
for j=1:N
    wave1 = [wave1 sound1(j,:)];
    wave2 = [wave2 sound2(j,:)];
    wave3 = [wave3 sound3(j,:)];
    wave4 = [wave4 sound4(j,:)];
    wave5 = [wave5 sound5(j,:)];
end

% save the sound to disk for subsequent playing 
audiowrite ('lab1.wav',wave1',fs); 
audiowrite ('lab2.wav',wave2',fs);
audiowrite ('lab3.wav',wave3',fs);
audiowrite ('lab4.wav',wave4',fs);
audiowrite ('lab5.wav',wave5',fs);

% repeat the above progress to produce 4 more sounds with more details

% Part G electrodes stimulu
% 4 is the number of electrodes
Etrodes=4; 
stimulus=zeros(N,Etrodes);
% N is the number of epochs
for i=1:N
    for j=1:5
        if   50<=ffreq(i,j) && ffreq(i,j)<750
            stimulus(i, 1) = amplitude(i);
        elseif 500<=ffreq(i,j) && ffreq(i,j)<1700
            stimulus(i,2)=amplitude(i);
        elseif 1500<=ffreq(i,j) && ffreq(i,j)<4000
            stimulus(i,3)=amplitude(i);
        elseif 3500<=ffreq(i,j) && ffreq(i,j)<11000
            stimulus(i,4)=amplitude(i);
        end
    end
end

% produce the eletrodogram plot
imagesc([1,N], [1:Etrodes], (stimulus'));
colorbar;

