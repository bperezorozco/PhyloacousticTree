%load mtlb;
tmp = get_random_TIER_file( filenames, '../../' )
[mtlb Fs] = audioread( tmp );
segmentlen= 100;
noverlap = 90;
NFFT = 128;

%spectrogram(mtlb, segmentlen, noverlap, NFFT, Fs, 'yaxis');

% dt = 1/Fs;
% I0 = round(0.1/dt);
% Iend = round(0.25/dt);
% 
% x = mtlb(I0:Iend);
x1 = mtlb.*hamming(length(mtlb));
x1 = filter(1, [1 0.63], x1);
x1 = mean_normalise(x1);
spectrogram(x1, segmentlen, noverlap, NFFT, Fs, 'yaxis');
colormap bone;

A = lpc(x1, 12);
rts = roots(A);
[h,f]=freqz(1,A,512,Fs);
figure;
plot(f,20*log10(abs(h)+eps));
legend('LP Filter');
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));

[frqs, indices] = sort(angz.*(Fs/(2*pi)), 'descend');
bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
res = [frqs bw]
nn = 1;
for kk = 1:length(frqs)
    if (frqs(kk) > 500 && bw(kk) <400)
        forms(nn) = frqs(kk);
        band(nn) = bw(kk);
    end
    nn = nn+1;
end