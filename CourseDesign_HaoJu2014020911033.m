close all;

PulseDuration=0.5;% 0.5s per code element
fs=1000;%Sampling Frequency: 1000 points/s
fc=50;%Carrier Frequency. Due to Nyquist Theroem,fc<500Hz
SNR=40;%SNR in dB
numCodeElement=PulseDuration*fs;%Number of samples within a pulse
N=8;%Length of Random Sequence
K=1;%Coefficient of the matched filter

TimeDuration=PulseDuration*N;%Total Time Duration
nofPulseDuration=1:PulseDuration*fs;%n sequence of a pulse
nofTimeDuration=1:TimeDuration*fs;%n sequence of the whole plot

% Random Sequence Generation: a random sequence with the length of N
Gate=0.5;%threshold of the sequence. In this Sequence P(0)=Gate, P(1)=1-Gate)
OriginalSignal=unifrnd(0,1,1,N);
Sequence=ones(1,N);
for n=1:N
    if OriginalSignal(n)>=Gate
        Sequence(n)=1;
    else
        Sequence(n)=0;
    end
end
SignalSource=rectpulse(Sequence,numCodeElement);%%Rectangular pulse shaping, with numCodeElement(PulseDuration*fs) samples per pulse
plot(nofTimeDuration/fs,SignalSource);
title('Original Sequence');
xlabel('time/s');

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Linear Frequency Modulation
a=1;%Amplitude of LFM Pulse
B=0.3;%Bandwidth of the LFM pulse
miu=B/(2*length(nofPulseDuration));
lfm=a*exp(j*(2*pi*fc*nofPulseDuration+2*pi*miu*nofPulseDuration.*nofPulseDuration));% LFM of a pulse
lfmwithNoise=awgn(lfm,SNR);%with White Gaussian Noise added
figure,plot(nofPulseDuration/fs,real(lfmwithNoise));%a Linear Frequency Modulation(LFM) pulse in time
title('LFM Pulse in Time');
xlabel('time/s');ylabel('Amplitude');
LFM=fftshift(fft(lfmwithNoise));%LFM pulse in Frequency
figure,plot(nofPulseDuration/TimeDuration,abs(LFM));
title('LFM Pulse in Frequency');
xlabel('freq/Hz');ylabel('Amplitude');

% Matched filter of LFM
nTimeDelay=1;
nt0=nTimeDelay+length(lfm);
nofMatchedFilter=1:length(lfm);
nofMatchedFilter=nt0-nofPulseDuration;
hLFM=K*lfm(nofMatchedFilter);%reversed part of the matched filter
hLFM=[zeros(1,nTimeDelay) hLFM];%matched filter with a beginning of a delay of TimeDelay/fs
figure,plot(linspace(1,PulseDuration,length(hLFM)),real(fliplr((hLFM))));
title('Matched Filter of LFM Pulse in Time');
xlabel('time/s');ylabel('Amplitude');
HLFM=fftshift(fft(hLFM));
figure,plot((1:length(HLFM))/PulseDuration,abs(HLFM));%matched filter in frequency 
% n transforming into w (or f?): frequency=n*(2pi/N)*(fs/2pi), N=TimeDuration*fs
title('Matched Filter of LFM Pulse in Frequency');
xlabel('freq/Hz');ylabel('Amplitude');

%Signal output of one pulse, after matched filtering
yLFM=conv(hLFM,lfmwithNoise);
tofyLFM=1:length(yLFM);
figure,plot(tofyLFM/fs,real(yLFM));%Signal output in time
title('LFM Pulse in Time, After Matched Filtering');
xlabel('time/s');ylabel('Amplitude');

%Linear phase modulation of the random sequence given above
SignalmodulatedbyLFM=zeros(1,TimeDuration*fs);
    for n1=0:(N-1)
        for n2=1:(PulseDuration*fs)
                SignalmodulatedbyLFM(PulseDuration*fs*n1+n2)=SignalSource(PulseDuration*fs*n1+n2)*lfm(n2);
        end
    end
SignalmodulatedbyLFM=awgn(SignalmodulatedbyLFM,SNR);   
figure,plot(nofTimeDuration/fs,real(SignalmodulatedbyLFM));
title('LFM Sequence in Time');
xlabel('time/s');

%LFM modulated signal passing through its corresponding matched filter
ySignalmodulatedbyLFM=conv(SignalmodulatedbyLFM,hLFM);
nofySignalmodulatedbyLFM=1:length(ySignalmodulatedbyLFM);
figure,plot(nofySignalmodulatedbyLFM/fs,real(ySignalmodulatedbyLFM));
title('Output of LFM Sequence after Matched Filtering');
xlabel('time/s');

%Recovery
ThresholdLFM=abs((0.5*a+a/(2*SNR)*log((1-Gate)/Gate))*max(ySignalmodulatedbyLFM))
tofRecovery=1:N;
LFMRecovery=zeros(1,N);
for Samples=1:N
    LFMRecovery(Samples)=ySignalmodulatedbyLFM(PulseDuration*fs*Samples);
    if LFMRecovery(Samples)>ThresholdLFM
        LFMRecovery(Samples)=1;
    else
        LFMRecovery(Samples)=0;
    end
end
LFMRecovered=rectpulse(LFMRecovery,numCodeElement);
figure,plot(nofTimeDuration/fs,LFMRecovered)
title('The LFM Sequence Recovered')
xlabel('time/s')

sumoflfm=0;
for n=1:length(lfm)
    sumoflfm=sumoflfm+real(lfm(n))*real(lfm(n));
end
N0=sumoflfm/(10.^(SNR/10));
SNRLFM=2*0.5*a*a*PulseDuration*fs/N0;
PeLFM=0.5*erfc(sqrt(SNRLFM/4))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ASK

a=1;%Amplitude of ASK=1;
carrier=a*cos(2*pi*fc*nofPulseDuration/fs);% the carrier signal of ASK
carrierwithNoise=awgn(carrier,SNR);%with White Gaussian Noise added
figure,plot(nofPulseDuration/fs,real(carrierwithNoise));%an ASK pulse in time
title('ASK Pulse in Time');
xlabel('time/s');
ASKPulse=fftshift(fft(carrierwithNoise));%LFM pulse in Frequency
figure,plot(nofPulseDuration/TimeDuration,abs(ASKPulse));
title('ASK Pulse in Frequency');
xlabel('freq/Hz');

%%%Matched filter 
nTimeDelay=1;
nt0=nTimeDelay+length(carrier);
nofMatchedFilterASK=nt0-nofPulseDuration;
hASK=K*carrier(nofMatchedFilterASK);%reversed part of the matched filter
hASK=[zeros(1,nTimeDelay) hASK];%matched filter with a beginning of a delay of TimeDelay/fs
figure,plot(linspace(1,PulseDuration,length(hASK)),real(hASK));%matched filter in time
title('Matched Filter of ASK Pulse in Time');
xlabel('time/s');
HASK=fftshift(fft(hASK));
figure,plot((1:length(HASK))/PulseDuration,abs(HASK));%matched filter in frequency 
%n transforming into w (or f?): frequency=n*(2pi/N)*(fs/2pi), N=TimeDuration*fs
title('Matched Filter of ASK Pulse in Frequency');
xlabel('freq/Hz');ylabel('Amplitude');

%Signal output of one pulse, after matched filtering
yASK=conv(hASK,carrierwithNoise);
tofyASK=1:length(yASK);
figure,plot(tofyASK/fs,real(yASK));%Signal output in time
title('ASK Pulse in Time, After Matched Filtering');
xlabel('time/s');

SignalModulatedbyASK=zeros(1,TimeDuration*fs);
    for n1=0:(N-1)
        for n2=1:(PulseDuration*fs)
                SignalModulatedbyASK(PulseDuration*fs*n1+n2)=SignalSource(PulseDuration*fs*n1+n2)*carrier(n2);
        end
    end
SignalModulatedbyASK=awgn(SignalModulatedbyASK,SNR);   
figure,plot(nofTimeDuration/fs,real(SignalModulatedbyASK));
title('ASK Sequence in Time');
xlabel('time/s');

SIGNALASK=fftshift(fft(SignalModulatedbyASK));
figure,plot(nofTimeDuration/TimeDuration,abs(SIGNALASK));
title('ASK Sequence in Frequency');
xlabel('Frequecy/Hz');

%%%ASK modulated signal passing through its corresponding matched filter
ySignalmodulatedbyASK=conv(SignalModulatedbyASK,hASK);
nofySignalmodulatedbyASK=1:length(ySignalmodulatedbyASK);
figure,plot(nofySignalmodulatedbyASK/fs,real(ySignalmodulatedbyASK));
title('Output of ASK Sequence after Matched Filtering');
xlabel('time/s');

%%%Recovery of ASK Signal
ThresholdASK=(0.5*a+a/(2*SNR)*log((1-Gate)/Gate))*max(ySignalmodulatedbyASK)
tofRecovery=1:N;
ASKRecovery=zeros(1,N);
for Samples=1:N
    ASKRecovery(Samples)=ySignalmodulatedbyASK(PulseDuration*fs*Samples);
    if ASKRecovery(Samples)>ThresholdASK
        ASKRecovery(Samples)=1;
    else
        ASKRecovery(Samples)=0;
    end
end
ASKRecovered=rectpulse(ASKRecovery,numCodeElement);
figure,plot(nofTimeDuration/fs,ASKRecovered)
title('The ASK Sequence Recovered')
xlabel('time/s')
% 
% sumofask=0;
% for n=1:length(lfm)
%     sumofask=sumofask+real(ask(n))*real(ask(n));
% end
% N0=sumofask/(10.^(SNR/10));
PeASK=0.5*erfc(sqrt(SNR/4))