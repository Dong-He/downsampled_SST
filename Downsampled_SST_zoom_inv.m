function [ Re_Sig ] = Downsampled_SST_zoom_inv( Tx,fs,freqr,win,hop,nfft,N,zoom,C)
%   inverse transform SST with subdivision
%   Note: the input parameters should be consistent with the forward
%   transform
%   Tx ： STFT or SST time-frequency matrix
%   win ： window function
%   hop: hop size of window (downsampling factor in time);
%   nfft: number of points for FFT
%   N：length of the input signal
%   C: the constant for amplitude correction

% this file is modified from 'ISTFT.m' originally by 
% Ilker Bayram, Istanbul Technical University,
% Modified by Zhaozhibin, Date: 2017.09, Email:zhaozhibin@stu.xjtu.edu.cn
% Modified by HeDong, Date: 2019.04, Email:hedong@stu.xjtu.edu.cn


Tx(isnan(Tx))=0;
L = length(win);

Nc= size(Tx,2);                   % get size

df=freqr(2)-freqr(1);
freq_0=(freqr(1)-df:-df:0);
Tx=[zeros(length(freq_0),Nc);Tx];

Tx=real(ifft(Tx,nfft*zoom));

Tx = Tx(1:L , :);
Tx = Tx.';
Re_Sig = zeros((Nc-1)*hop+L , 1);

Half = floor(L/2);

for n = 1 : L
    index = (n : hop :(Nc-1)*hop+n);
    Re_Sig(index) = Re_Sig(index) + win(n) .* Tx(: , n);
end

Re_Sig = Re_Sig(Half+1:Half+N)*nfft*zoom/fs*hop/C;   % 注意这里要乘以2
% Re_Sig = Re_Sig(Half+1:Half+N)*nfft*zoom/fs*hop;   % 注意这里要乘以2

end

