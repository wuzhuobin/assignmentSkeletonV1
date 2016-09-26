function rrcisim(nch,filename,bpflow,recttype,lpftype,outdir)
% Usage: rrcisim (nch,filename,bpflow,recttype,lpftype,outdir)
%
% nch -- total number of channels
% filename - Input filename (you need a directory called "audio_in" for input audio files)
% bpflow - cutoff of lowest BPF (90 or 350Hz)
% recttype - rectification ('f' - full, 'h' - half)
% lpftype - lowpass filter type ('mvav' - 16pt moving average, 'b400' - 400Hz 2nd order butterworth
% outdir - directory and filename ex. 'real_instruments/piano' (you need a directory called "audio_out")
% -- calls function estfilt.m for filter parameters --
%
% adapted from Loizou (JASA 1999)
%
% Rebecca Reich
% January 2002
global filterA filterB center Srate
%========== open the input file ========== %
% ====== USING WAVREAD (so no need to scale output) ========
[x,Srate,bpsa] = wavread(['audio_in/', filename]);
[row,col] = size(x);
if col==2d
x = sum(x,2);
 disp('summing to mono');
end
% now resample to 12971Hz
if Srate==44100
 x = resample(x,5,17);
 Srate = Srate*5/17;
elseif Srate == 22050
 x = resample(x,10,17);
 Srate = Srate*10/17;
elseif Srate == 48000
 x = resample(x,459375,1700000);
 Srate = Srate*459375/1700000;
end
disp('downsample to 12971Hz');
% ===============================================
n_samples=length(x);
nChannels=nch;
% ========== remove any DC bias ==========
x = x - mean(x);
% %--Preemphasize first (LOIZOU method) ---------
% bp = exp(-1200*2*pi/Srate);
% ap = exp(-3000*2*pi/Srate);
% x = filter([1 -bp],[1 -ap],x); %using freqz: cuts off high freqs?
% ========== get bandpass filter coefficients ==========
if isempty(filterA)
 estfilt(nChannels,bpflow);
 fprintf('\n Getting filter coefficients for %d-channel processor..\n',nch);
end
% ========== calculate lowpass filter coefficients ==========
if lpftype == 'mvav'
mvavgord = 16; % sample length of moving average filter (LPF smoothing filter)
blow = 1/(mvavgord) * ones(1,mvavgord);
alow = 1;
elseif lpftype == 'b400'
 [blow,alow]=butter(2,400/(Srate/2)); % b400 has 2nd order filter for signals done at 44.1/22.05kHz
 %[blow,alow] = butter(6,400/(Srate/2));
else
 [blow,alow]=butter(2,40/(Srate/2)); % an attempt at a lower cutoff frequency
end
% ========== filter input with BPFs, lowpass filter the rectified output ==========
ylpf=zeros(nChannels,n_samples);
ybpf=zeros(nChannels,n_samples);
for i=1:nChannels
 ybpf(i,:)=filter(filterB(i,:),filterA(i,:),x)'; % bandpass
 if recttype == 'f'
 yrect(i,:) = abs(ybpf(i,:)); % full-wave rectify
 ylpf(i,:)=filter(blow,alow,yrect(i,:)); % lpf
 else % half-wave
 yrect(i,:) = ybpf(i,:);
 yrect(i,yrect(i,:)<0) = 0; % set neg. values to zero
 ylpf(i,:) = filter(blow,alow,yrect(i,:)); % lpf
 end
end
%ysum = sum(ylpf,1); %for playing purposes
disp('done filtering');
% HIGHPASS (LOIZOU method) ----------------
% [bhi,ahi] = butter(4,20/(Srate/2));
% for i = 1:nChannels
% ylpf(i,:) = filter(bhi,ahi,ylpf(i,:));
% end
% ------------------------
% ========== outputs ==========
freq=center/Srate; %normalized center frequencies
for i=1:nChannels % ------ generate sinewaves ---------
 ycos(i,:)=cos(2*pi*freq(i)*[0:n_samples-1]);
 yout_off(i,:)=ylpf(i,:).*ycos(i,:); %modulate sinewaves with envelopes
 yout(i,:) = yout_off(i,:);
end
youtsum = sum(yout,1);
youtsum = youtsum/(max(abs(youtsum)));
if max(abs(youtsum))>1, fprintf('Warning! Overflow in file %s\n',filename); end;
%--- save output to a file ----
wavname = ['_ch',num2str(nch),'_BPF',num2str(bpflow),'_RECT',recttype,'_LPF',num2str(lpftype),'.wav'];
name = ['audio_out/',outdir,wavname]
wavwrite(youtsum,Srate,bpsa,['audio_out/',outdir,wavname]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME AND FREQUENCY VECTORS FOR PLOTTING
% tvec = 0:1/Srate:(n_samples-1)/Srate;
% fvec = [0:n_samples-1]./n_samples*Srate;