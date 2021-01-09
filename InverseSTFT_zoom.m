function [ Re_Sig ] = InverseSTFT_zoom( X , Win ,fs, freqr,Hop_Size ,N,Nfft)
% Inverse STFT
% this file is modified from 'STFT.m' originally by 
% Ilker Bayram, Istanbul Technical University,
% Input:
%     X : The short time fourier transforming coefficients
%     Win : The length of the window
%     Hop_Size : the size of the jump
%     N : the length of the signal
% Output:
%     Re_Sig : the reconstructed signal
% Attention : the length of FFT is equal to the length of the Windows
% Modified by Zhaozhibin
% Date: 2015.09
% Email:zhaozhibin@stu.xjtu.edu.cn
X(isnan(X))=0;
L = length(Win);
[~, Nc] = size(X);                   % get size

df=freqr(2)-freqr(1);
freq_0=(freqr(1)-df:-df:0);
X=[zeros(length(freq_0),Nc);X];

c = real(ifft(X,Nfft));
c = c(1:L , :);
c = c.';
Re_Sig = zeros((Nc-1)*Hop_Size+L , 1);

Half = floor(L/2);
for n = 1 : L
    index = (n : Hop_Size :(Nc-1)*Hop_Size+n);
    Re_Sig(index) = Re_Sig(index) + Win(n) .* c(: , n);
end
% Re_Sig = Re_Sig(Half+1:Half+N) * sqrt(Nfft) *  Hop_Size;
Re_Sig = Re_Sig(Half+1:Half+N) * Nfft/fs *  Hop_Size;

end

