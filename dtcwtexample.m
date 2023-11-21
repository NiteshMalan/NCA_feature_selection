clc;
clear all;
close all;
load('c4.mat');
load('c3.mat');
load('eyecentral.mat');
c3=c3(106251:108250);
c4=c4(106251:108250);
eyecentral=eyecentral(106251:108250);
fs=250;
N=length(c3);

M=round(fs/3.7); %ssa window length
t=(0:1:N-1)/fs;

%calculating input SNR
c3=c3-mean(c3);
% c3=c3/std(c3,1);
eyecentral=(eyecentral-mean(eyecentral));
% eyecentral=eyecentral/std(eyecentral,1)
avg_eeg=sum((abs(c3)).^2);
avg_eog=sum((abs(eyecentral)).^2);
X = c3; %+eyecentral;
X=X-mean(X);
[peaksnr, snr] = psnr(X, c3);
fprintf('\n The Input Peak-SNR value is %0.4f', peaksnr);
fprintf('\n The  Input SNR value is %0.4f \n', snr);
snr_input= 10*log(avg_eeg/avg_eog);
figure(1), clf;
% subplot(311);
plot(t,c3,'red');
set(gca,'fontsize', 24)
axis tight
title('EEG signal','fontsize',24);
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Time(in sec)','fontsize',24)
% ylim([-50 250])
% subplot(312);
% plot(t,eyecentral,'blue');
% title('EOG Signal','fontsize',12);
% ylabel('Ampluted (in uV)','fontsize',12)
% xlabel('Time(in sec)','fontsize',12)
% ylim([-50 250])
% subplot(313);
% plot(t,X,'green');
% title('Contaminated EEG signal','fontsize',12);
% ylabel('Ampluted (in uV)','fontsize',12)
% xlabel('Time(in sec)','fontsize',12)
% % ylim([-50 250])
% axis tight
[Faf, Fsf] = FSfarras; % 1st stage anal. & synth. filters
[af, sf] = dualfilt1;
wt = dddtree('cplxdt',X,4,Faf,af);

outputindices = {[1,1];[1,2];[2 1]; [2 2]; [3 1]; [3 2];[4,1];[4,2];[5,1];[5,2]};
out = dddtreecfs('e',wt,'ind',outputindices);

coef1=out{1,1}+i.*out{1,2};
coef2=out{1,3}+i.*out{1,4};
coef3=out{1,5}+i.*out{1,6};
coef4=out{1,7}+i.*out{1,8};
coef5=out{1,9}+i.*out{1,10};

artifact = cell2mat(dddtreecfs('r',wt,'scale',{5}));

output_sig = dddtreecfs('r',wt);
figure (2);

plot(artifact,'r');
title('Contaminated EEG and DTCWT decomposed at sacle 5')
axis tight;

output=X-artifact;

dtcwt_avgp_output=sum(abs(output).^2)/N;
dtcwt_MSE=(sum((c3-output).^2))/N;
dtcwt_snr_output = 10*log(dtcwt_avgp_output/dtcwt_MSE);
dtcwt_rmsse=(sum((c3-output).^2))/sum((c3).^2);
[peaksnr1, snr1] = psnr(output, c3);
fprintf('\n The DTCWT Peak-SNR value is %0.4f', peaksnr1);
fprintf('\n The  DTCWT SNR value is %0.4f \n', snr1);
fprintf('\n The  DTCWT RRMSE value is %0.4f \n', dtcwt_rmsse);
noise_var= (median(abs(artifact-median(abs(artifact))))/0.6745);

T = noise_var*(sqrt(2*log(N)));
%xWT = perform_thresholding(wt,T,'hard');
%soft threhold
soft_thres1 = wthresh(coef1,'s',T);
soft_thres2 = wthresh(coef2,'s',T);
soft_thres3 = wthresh(coef3,'s',T);
soft_thres4 = wthresh(coef4,'s',T);
soft_thres5 = wthresh(coef5,'s',T);

st=wt;
st.cfs{1,5}(:,:,1)=soft_thres5.*real(coef5);
st.cfs{1,5}(:,:,2)=soft_thres5.*imag(coef5);
st.cfs{1,4}(:,:,1)=soft_thres4.*real(coef4);
st.cfs{1,4}(:,:,2)=soft_thres4.*imag(coef4);
st.cfs{1,3}(:,:,1)=soft_thres3.*real(coef3);
st.cfs{1,3}(:,:,2)=soft_thres3.*imag(coef3);
st.cfs{1,2}(:,:,1)=soft_thres2.*real(coef2);
st.cfs{1,2}(:,:,2)=soft_thres2.*imag(coef2);
st.cfs{1,1}(:,:,1)=soft_thres1.*real(coef1);
st.cfs{1,1}(:,:,2)=soft_thres1.*imag(coef1);

soft_RC=(dddtreecfs('r',st));

figure(41);
plot(t,soft_RC)
title('soft threshold Extracted Output Noise')

soft_out_eeg=X-soft_RC;
figure(61)
plot(t,soft_out_eeg) 
title('soft threshold Corrected Output EEG')

k_omega5=zeros(1,125);
b=T/10
for i=1:125;
k_omega5(i)=1/(1+exp(-(abs(coef5(i))-T)/b));

R5=[k_omega5];
out_noise5(i)=[1-k_omega5(i)];
% coe=(R5).*0+(out_noise5).*1;
end
k_omega4=zeros(1,125);
for j=1:125;
k_omega4(j)=1/(1+exp(-(abs(coef4(j))-T)/b));
R4=[k_omega4];
out_noise4(j)=[1-k_omega4(j)];
end

k_omega3=zeros(1,250);
for k=1:250;
k_omega3(k)=1/(1+exp(-(abs(coef3(k))-T)/b));
R3=[k_omega3];
out_noise3(k)=[1-k_omega3(k)]; 
end

k_omega2=zeros(1,500);
for l=1:500;
k_omega2(l)=1/(1+exp(-(abs(coef2(l))-T)/b));
R2=[k_omega2];
out_noise2(l)=[1-k_omega2(l)];
end

k_omega1=zeros(1,1000);
for m=1:1000;
k_omega1(m)=1/(1+exp(-(abs(coef1(m))-T)/b));
R1=[k_omega1];
out_noise1(m)=[1-k_omega1(m)];
end

figure(3);
subplot(511)
plot(R5)
subplot(512)
plot(R4)
subplot(513)
plot(R3)
subplot(514)
plot(R2)
subplot(515)
plot(R1)


rt=wt;
rt=wt;
rt.cfs{1,5}(:,:,1)=R5.*real(coef5);
rt.cfs{1,5}(:,:,2)=R5.*imag(coef5);
rt.cfs{1,4}(:,:,1)=R4.*real(coef4);
rt.cfs{1,4}(:,:,2)=R4.*imag(coef4);
rt.cfs{1,3}(:,:,1)=0;
rt.cfs{1,3}(:,:,2)=0;
rt.cfs{1,2}(:,:,1)=0;
rt.cfs{1,2}(:,:,2)=0;
rt.cfs{1,1}(:,:,1)=0;
rt.cfs{1,1}(:,:,2)=0;


RC=(dddtreecfs('r',rt));

figure(4);
subplot(211)
plot(t,RC)
title('QAWS Extracted Output Noise')
subplot(212)
plot(t,eyecentral)

figure(5);
subplot(411)
plot(abs(coef4))
subplot(412)
plot(abs(coef3))
subplot(413)
plot(abs(coef2))
subplot(414)
plot(abs(coef1))

final=X-RC;
figure(6)
% subplot(211)
hold on
plot(t,final) 
set(gca,'fontsize', 24)
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Time(in sec)','fontsize',24)
axis tight
% ylim([-50 250])
title('DTCWT-QWS Corrected Output EEG','fontsize',24)
plot(t,c3,'red')
set(gca,'fontsize', 24)
axis tight
legend('DTCWT-QWS Corrected Output EEG','Real EEG')
% subplot(212)
% plot(t,RC)
% set(gca,'fontsize', 18)
% ylabel('Ampluted (in uV)','fontsize',18)
% xlabel('Time(in sec)','fontsize',18)
% ylim([-70 250])
% title('Extracted EOG signal using DTCWT-QWS','fontsize',24)
figure(7), clf, 
subplot(211)
plot(t,final)
ylim([-50 100])
subplot(212)
plot(t,c3,'r');
ylim([-50 100])
avgp_output=sum(abs(final).^2)/N;
MSE=(sum((c3-final).^2))/N;
snr_output = 10*log(avgp_output/MSE);
rmsse=(sum((c3-final).^2))/sum((c3).^2);
[peaksnr2, snr2] = psnr(final, c3);
fprintf('\n The DTCWT-QWS Peak-SNR value is %0.4f', peaksnr2);
fprintf('\n The  DTCWT-QWS SNR value is %0.4f \n', snr2);
fprintf('\n The  DTCWT-QWS RRMSE value is %0.4f \n', rmsse);

figure(71)
fbins=(0:1/N:1-1/N)*fs;
subplot(122)
plot(fbins,abs(fft(RC)))
set(gca,'fontsize', 24)
xlim([0 30])
title('FFT of the DTCWT-QWS extracted EOG','fontsize',24)
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Freq (in Hz)','fontsize',24)
subplot(121)
plot(fbins,abs(fft(eyecentral)))
set(gca,'fontsize', 24)
xlim([0 30])
title('FFT of the added EOG signal','fontsize',24)
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Freq (in Hz)','fontsize',24)

figure(72)
fbins=(0:1/N:1-1/N)*fs;
subplot(122)
plot(fbins,abs(fft(final)))
set(gca,'fontsize', 18)
xlim([0 30])
title('FFT of the DTCWT-QWS extracted EEG','fontsize',24)
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Freq (in Hz)','fontsize',24)
subplot(121)
plot(fbins,abs(fft(c3)))
set(gca,'fontsize', 18)
xlim([0 30])
title('FFT of the added EEG signal','fontsize',24)
ylabel('Ampluted (in uV)','fontsize',24)
xlabel('Freq (in Hz)','fontsize',24)

[c,l] = wavedec(X,7,'db3');
[cd1,cd2,cd3,cd4,cd5,cd6,cd7] = detcoef(c,l,[1 2 3 4 5 6 7]);

k_omega1dwt=zeros(1,2031);
for g=1:2031;
k_omega1dwt(g)=1/(1+exp(-(abs(c(g))-T)/b));
Rdwt=[k_omega1dwt];
out_noisedwt(g)=[1-k_omega1dwt(g)];
end
dwtnoise=c.*Rdwt;
% dwteeg=c-c.*Rdwt;
a1 = waverec(c.*Rdwt,l,'db3');
%a2 = waverec(dwteeg,l,'db3');
dwtout=X-a1;
figure(200)
subplot(211)
plot(t,a1)
subplot(212)
plot(t,dwtout)
rmsse4=(sum((c3-dwtout).^2))/sum((c3).^2);

[peaksnr4, snr4] = psnr(dwtout, c3);
fprintf('\n The DWT-QWS threshold Peak-SNR value is %0.4f', peaksnr4);
fprintf('\n The  DWT-QWS threshold value is %0.4f \n', snr4);
fprintf('\n The  DWT-QWS RRMSE value is %0.4f \n', rmsse4);

% cd = wthresh(c,'s',T);
% a0 = waverec(cd,l,'db3');

% figure(89)
% clf;
% dwt_hard=perform_thresholding(c,T,'hard')
% dwt_soft=perform_thresholding(c,T,'soft');
% dwt_semisoft1=perform_thresholding(c,[T 2*T],'semisoft');
dwt_semisoft2=perform_thresholding(c,[T 4*T],'semisoft');
a0 = waverec(dwt_semisoft2,l,'db3');
figure(101)
plot(t,a0)

final1=X-a0;
figure(106)
plot(t,final1) 
title('DWT Corrected Output EEG')
figure(107), clf, 
subplot(211)
plot(t,final1)
ylim([-50 100])
subplot(212)
plot(t,c3,'r');
ylim([-50 100])
avgp_output1=sum(abs(final1).^2)/N;
MSE1=(sum((c3-final1).^2))/N;
snr_output1 = 10*log(avgp_output1/MSE1);
rmsse1=(sum((c3-final1).^2))/sum((c3).^2);

[peaksnr3, snr3] = psnr(final1, c3);
fprintf('\n The DWT-Soft threshold Peak-SNR value is %0.4f', peaksnr3);
fprintf('\n The  DWT-Soft threshold value is %0.4f \n', snr3);
fprintf('\n The  DWT-soft RRMSE value is %0.4f \n', rmsse1);

%% morlet wavelet
srate = fs; % in hz
time  = -0.5:1/srate:0.5; % best practice is to have time=0 at the center of the wavelet
frex  = 10; % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp( 1i*2*pi*frex.*time );

% create Gaussian window
s = 7 / (2*pi*frex); % this is the standard deviation of the gaussian
gaus_win  = exp( (-time.^2) ./ (2*s^2) );


% now create Morlet wavelet
cmw = sine_wave.*gaus_win;
kernel=real(cmw);

nkernel =length(kernel);
nconv = N+nkernel-1;

Xcmw= fft(cmw,nconv);
Xcmw = Xcmw./max(Xcmw);
Xsignal = fft(final,nconv);
conv_res =Xsignal'.*Xcmw';

hz = linspace(0,srate/2,floor(nconv/2)+1);

figure(8), clf

% plot power spectrum of data
subplot(311)
plot(hz,2*abs(Xsignal(1:length(hz))/length(output)))

% plot power spectrum of wavelet
subplot(312)
plot(hz,abs(Xcmw(1:length(hz))))

% plot power spectrum of convolution result
subplot(313)
plot(hz,2*abs(conv_res(1:length(hz))/length(output)))

% now timedomain
length(conv_res);
% cut 1/2 of the length of the wavelet from the beginning and from the end
half_wav = floor( length(cmw)/2 )+1;

% take inverse Fourier transform
conv_res_timedomain = ifft(conv_res);

conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav);

figure(9), clf
subplot(211)
plot(t,output,'k')
subplot(212)
plot(t,real(conv_res_timedomain),'r')
legend({'EEG data';'convolution-filtered data'})

%% fft of noise and extracted noise
fbins= (0:1/N:1-1/N)*fs;
eog_fft= fft(eyecentral);
abs_eog_fft = abs(eog_fft);
fft_RC= fft(RC);
abs_fft_RC= abs(fft_RC);
figure(10)
subplot(211)
plot(fbins,abs_eog_fft)
subplot(212)
plot(fbins,abs_fft_RC)
figure(11)
spectrogram(final,[],[],[],fs,'yaxis')

COEFS = cwt(final,1:64,'cgau4');

% Compute and plot the scalogram (image option)
figure(12);
SC = wscalogram('image',COEFS);
figure(13), plot(t, SC);
cfreq = centfrq('sym4');

