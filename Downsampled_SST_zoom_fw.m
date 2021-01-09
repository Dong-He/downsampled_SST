function [s,Tx,time,freqr_stft,freqr_sst,win]=Downsampled_SST_zoom_fw(sig,fs,sigma,hop,f1,f2,zoom,nfft)
% downsampled synchrosqueezing transform, with reconstruction property
% sig: input signal£»fs£ºsampling frequency£»sigma£ºGaussian window parameter£»
% hop: hop size of window (downsampling factor in time); nfft: number of points for FFT
% [f1,f2]: the freq range; zoom: the subdivision factor

% this file is modified from 'ISTFT.m' originally by 
% Ilker Bayram, Istanbul Technical University,
% Modified by Zhaozhibin, Date: 2017.09, Email:zhaozhibin@stu.xjtu.edu.cn
% Modified by HeDong, Date: 2019.04, Email:hedong@stu.xjtu.edu.cn

if nargin < 8
    nfft=length(sig);
end

sig=sig(:);

idx1=1+round(f1*nfft/fs); %f1 index
idx2=1+round(f2*nfft/fs); %f2 index
freqr=(0:nfft-1)*fs/nfft;
% Gaussian window function
K = 0.005;
Half = round(sqrt(-2*log(K))*sigma*fs); 
sig = [zeros(Half,1); sig; zeros(Half-hop+1,1)];  % padded signal

% window function
ix     = ((-Half):Half)';
t_win = ix/fs;
win = (pi*sigma^2)^(-0.25).* exp(-(t_win/sigma).^2/2);
% the first derivative of Gaussian window function
dwin = (pi*sigma^2)^(-0.25).* exp(-(t_win/sigma).^2/2).*(-t_win/sigma^2);
L=length(win);

% do STFT
ds=buffer(sig,L,L-hop, 'nodelay');
s = bsxfun(@times, win(:), ds);
s = fft(s,nfft)*2/nfft;
s = s(idx1:idx2,:);
ds = bsxfun(@times, dwin(:), ds);
ds = fft(ds,nfft)*2/nfft;
ds = ds(idx1:idx2,:);

gamma = 1e-8;       
s((abs(s)<gamma)) = NaN;
% make time and freq axis
[row,col]=size(s);
time=(0:hop:(col-1)*hop)/fs;
freqr_stft= freqr(idx1:idx2)';

CandidateIF =  real(1i *ds ./s / (2*pi) ) ;  
CandidateIF = bsxfun(@plus, freqr_stft(:), CandidateIF);

% multiply by a phase term for the following reassignment operation
s  = bsxfun(@times, exp(1i*2*pi*Half/fs*freqr_stft), s);
freqr_sst= linspace(freqr(idx1),freqr(idx2),zoom*(idx2-idx1)+1)';  % update freq axis
Tx = SynchroSqueezing(s,CandidateIF,freqr_sst);
% multiply back by a inverse phase term for the following reconstruction
Tx  = bsxfun(@times, exp(-1i*2*pi*Half/fs*freqr_sst), Tx);
s  = bsxfun(@times, exp(-1i*2*pi*Half/fs*freqr_stft), s);
end

% reassignment operation
function Tx = SynchroSqueezing(STFT,CandidateIF,freqr)
% STFT-based Synchrosqueezing 
%   input:
%       STFT: the TF representation by STFT
%       w: Candidate instantaneous frequency
%       freqr: the frequency associated with STFT
%   output:
%       Tx: `the synchrosqueezing result

% if nargin ~= 3
%     error('you should check the input parameter of SynSqu_STFT function');
% end
[STFT_rows,STFT_cols] = size(STFT);
Tx = zeros(length(freqr),STFT_cols);  
delta_f = (freqr(2)-freqr(1));  % STFT_rows is No. freq bins, STFT_cols is No. time bins
k=zeros(size(Tx));
for u=1:STFT_cols
   for fi=1:STFT_rows
        if (~isnan(CandidateIF(fi, u)) && (CandidateIF(fi,u)>0))
            % Find w_l nearest to w(a_i,b)         
            k(fi,u) = 1+round((CandidateIF(fi,u)-freqr(1))/delta_f);
            if ~isnan(k(fi,u)) && k(fi,u)>0 && k(fi,u) <= size(Tx,1)
                Tx(k(fi,u),u) = Tx(k(fi,u), u) + STFT(fi, u); 
            end
        end
    end 
end 

end