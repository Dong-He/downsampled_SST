function [ Re_Sig] = Downsampled_inv( Tx,fs,win,hop,nfft,N,C)
%   inverse transform for STFT and SST
%   Note: the input parameters should be consistent with the forward
%   transform
%   Tx £º STFT or SST time-frequency matrix
%   win £º window function
%   hop: hop size of window (downsampling factor in time);
%   nfft: number of points for FFT
%   N£ºlength of the input signal
%   C: the constant for amplitude correction

% this file is modified from 'ISTFT.m' originally by 
% Ilker Bayram, Istanbul Technical University,
% Modified by Zhaozhibin, Date: 2017.09, Email:zhaozhibin@stu.xjtu.edu.cn
% Modified by HeDong, Date: 2019.04, Email:hedong@stu.xjtu.edu.cn

Tx(isnan(Tx))=0;
L = length(win);

[~, Nc] = size(Tx);                   % get size

Tx=real(ifft(Tx,nfft));

Tx = Tx(1:L , :);
Tx=bsxfun(@times,Tx,win(:));

Re_Sig = zeros((Nc-1)*hop+L , 1);
% winSig = zeros((Nc-1)*hop+L , 1);
Half = floor(L/2);
% A=0;

for n = 1 : L
    index = (n : hop :(Nc-1)*hop+n);
%     Re_Sig(index) = Re_Sig(index) + win(n) .* Tx(: , n);
Re_Sig(index) = Re_Sig(index) + Tx(n,:)';
%     winSig(index) = winSig(index) + win(n).^2;

%     A(n)=sum(win(index).^2);
end
% figure
% plot(winSig/fs*hop)
% Re_Sig = Re_Sig(Half+1:Half+N) * sqrt(nfft) *  hop;
Re_Sig = Re_Sig(Half+1:Half+N)*nfft/fs*hop/C;
% Re_Sig = Re_Sig(Half+1:Half+N)*nfft/fs*hop/sqrt(2);

end

