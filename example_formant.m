load mtlb;
%[mtlb Fs] = audioread('../../dataset/Acrocephalus arundinaceus/AcrAru00002.wav');
segmentlen= 100;
noverlap = 90;
NFFT = 128;

%spectrogram(mtlb, segmentlen, noverlap, NFFT, Fs, 'yaxis');

dt = 1/Fs;
I0 = round(0.1/dt);
Iend = round(0.25/dt);

x = mtlb(I0:Iend);
x1 = x.*hamming(length(x));
x1 = filter(1, [1 0.63], x1);
spectrogram(x1, segmentlen, noverlap, NFFT, Fs, 'yaxis');

A = lpc(x1, 12)
rts = roots(A);

rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));

[frqs, indices] = sort(angz.*(Fs/(2*pi)), 'descend')
bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)))
%res = [frqs bw]
% nn = 1;
% for kk = 1:length(frqs)
%     forms(nn) = frqs(kk);
%     nn = nn+1;
%     %if (frqs(kk) > 90 && bw(kk) <400)
%         
%     %end
% end
% forms