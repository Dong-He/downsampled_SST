% 
% Copyright (C) 2019 He Dong <hedong@stu.xjtu.edu.cn>

close all;clear;clc;

fs=1024;
t=(0:1024*8-1)/fs;
N=length(t);
data(1,:)=(1-0.1*cos(0.25*pi*t)).*cos(100*pi*t);
t1=t(1:N/2-1);
t2=t(N/2:end);
data(2,:)=[cos(500*pi*t1-25*pi*t1.^2),cos(500*pi*t2-50*pi*t2.^2+25/6*pi*t2.^3+pi*4/3)];
data(3,:)=(1-0.2*cos(pi/8*t)).*cos(740*pi*t+50*sin(1*pi*t));
sig=sum(data);
randn('state',1);   
x=awgn(sig,5,'measured');   
sigma=0.03;

%% plot theoretical IF
f1=50*ones(size(t));
f2=[250-25*t1,250-50*t2+6.25*t2.^2];
f3=370+50*cos(0.75*pi*t)-50*cos(0.5*pi*t);

figure
plot(t,f1,'k')
hold on
plot(t,f2,'k')
hold on
plot(t,f3,'k')


%% downsampled SST

nfft=length(x);
hop=20;
[~,Tx,time,freqr,~]=Downsampled_SST_fw(x,fs,sigma,hop,nfft);
figure
imagesc(time,freqr,abs(Tx));shading interp;axis xy

%% selective reassignment: x1 component[0,80]Hz
[~,Tx,time,freqr,~]=Downsampled_SST_zoom_fw(x,fs,sigma,hop,0,80,1,nfft);
figure
imagesc(time,freqr,abs(Tx));shading interp;axis xy


%% reconstruction

nfft=length(x)/4;
hop=20;
data=data.';    

% extract TF ridges
[s,Tx,time,freqr,win]=Downsampled_SST_fw(x,fs,sigma,hop,nfft);
winlen=length(win);
C=sum(win((winlen-1)/2)*win)/fs;
nb=round(10/(freqr(2)-freqr(1)));
penval=0.001;
[IF_sst,idx_sst] = tfridge(abs(Tx),freqr,penval,'NumRidges',3,'NumFrequencyBins',nb);
[IF_stft,idx_stft] = tfridge(abs(s),freqr,penval,'NumRidges',3,'NumFrequencyBins',nb);
k=1;
for wd=20
    Tx1=Tx;
    s1=s;
    nbins=max(round(wd/(freqr(2)-freqr(1))),1);
    for i=1:3
        
        H_sst=zeros(size(Tx1));
        H_stft=zeros(size(s1));
        for j=1:size(Tx1,2)
            H_sst(idx_sst(j,i)-nbins:idx_sst(j,i)+nbins,j)=1;
            H_stft(idx_stft(j,i)-nbins:idx_stft(j,i)+nbins,j)=1;
        end
        Txf=Tx1.*H_sst;
        sf=s1.*H_stft;
        Re_Sig_sst(:,i,k) = Downsampled_inv(Txf,fs,win,hop,nfft,length(x),C);
        Re_Sig_stft(:,i,k) = Downsampled_inv(sf,fs,win,hop,nfft,length(x),1);
        s1=s1-sf;
        Tx1=Tx1-Txf;
    end
    
    % get SNR
    
    e_sst(:,:,k)=bsxfun(@minus,Re_Sig_sst(:,:,k),data);
    e_stft(:,:,k)=bsxfun(@minus,Re_Sig_stft(:,:,k),data);
    for i=1:3
        S_sst(k,i)=snr(data(:,i),e_sst(:,i,k));
        S_stft(k,i)=snr(data(:,i),e_stft(:,i,k));
    end
    
    k=k+1;
end
figure
imagesc(abs(s));axis xy;

figure
imagesc(abs(Tx));axis xy;

figure
x_sst=Re_Sig_sst(:,3,1);
x_stft=Re_Sig_stft(:,3,1);
plot(x_sst)
hold on
plot(x_stft)



