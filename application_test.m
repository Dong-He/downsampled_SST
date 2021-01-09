% 
% Copyright (C) 2019 He Dong <hedong@stu.xjtu.edu.cn>
% Application of downsampled SST on the large scale aero-engine data analysis
close all;clear;clc;

load('data3.mat')

sigma=0.02;

%% fast SST forward transform for large scale data

nfft=2^13;
L=1+2*round(sqrt(-2*log(0.005))*sigma*fs);
overlap=0.8;
hop=round((1-overlap)*L);
zoom=8;

f1=6900;
f2=9500;

[S,Tx,time,freqr_stft,freqr_sst,win]=Downsampled_SST_zoom_fw(x,fs,sigma,hop,f1,f2,zoom,nfft);

f_start=8500;
f_end=9500;

df=freqr_sst(2)-freqr_sst(1);
idx_f1=round((f_start-f1)/df);
idx_f2=round((f_end-f1)/df);

freqr_sst=freqr_sst(idx_f1:idx_f2);

t_start=100;
dt=time(2)-time(1);
idx_t1=round(t_start/dt);
% t2=end;


Tx=Tx(idx_f1:idx_f2,idx_t1:end); 
figure
imagesc(time(idx_t1:end),freqr_sst,abs(Tx));axis xy;shading interp;
figure
imagesc(time,freqr_stft,abs(S));axis xy;shading interp;
xlim([t_start,time(end)]);
ylim([f_start,f_end]);
wd=10;   % half bandwidth
penval=0.1;
nbins=round(wd/df);


%% TF ridge extraction
[IF_sst,idx_sst] = tfridge(abs(Tx),freqr_sst,penval,'NumRidges',1,'NumFrequencyBins',nbins);

idx_f1=round((f_start-f1)/(freqr_stft(2)-freqr_stft(1)));
idx_f2=round((f_end-f1)/(freqr_stft(2)-freqr_stft(1)));

freqr_stft=freqr_stft(idx_f1:idx_f2);

S=S(idx_f1:idx_f2,idx_t1:end);   

[IF_stft,idx_stft] = tfridge(abs(S),freqr_stft,penval,'NumRidges',1,'NumFrequencyBins',nbins);


%% TF recon for STFT
t=(0:length(x)-1)/fs;
wd=20;   % half bandwidth
nbins=round(wd/(freqr_stft(2)-freqr_stft(1)));
H_stft=zeros(size(S));

for j=1:size(S,2)
    H_stft(idx_stft(j)-nbins:idx_stft(j)+nbins,j)=1;
end
Sf=S.*H_stft;
N=round((t(end)-t_start)*fs);
Re_Sig_stft= InverseSTFT_zoom( Sf , win, fs,freqr_stft, hop,N,nfft);

t_recon=(0:N-1)/fs+t_start;

%% TF recon for SST
wd=10;   % half bandwidth
nbins=round(wd/(freqr_sst(2)-freqr_sst(1)));
H_sst=zeros(size(Tx));
for j=1:size(Tx,2)
    H_sst(idx_sst(j)-nbins:idx_sst(j)+nbins,j)=1;
end
Txf=Tx.*H_sst;
N=round((t(end)-t_start)*fs);

clear p S Sf Tx
C=sum(win((length(win)-1)/2)*win)/fs;
Re_Sig_sst=Downsampled_SST_zoom_inv( Txf,fs,freqr_sst,win,hop,nfft,N,zoom,C);

figure
plot(t_recon,Re_Sig_stft)
hold on
plot(t_recon,Re_Sig_sst)

%% compare recon RMSE of STFT and SST

RMSE=sqrt(mean((Re_Sig_sst-Re_Sig_stft).^2))

RMS_STFT=sqrt(mean(Re_Sig_stft.^2))



