function [ frqs ] = formants_from_mfcc( y, fs, P )
%SYNTHETIC_FROM_MFCC Summary of this function goes here
%   Detailed explanation goes here

fan = fs/2;
Man = freq2mel(fan);
NFFT = 2^nextpow2( length(y) );
f = mel2freq( Man * linspace(0, 1, NFFT/2 + 1) );
r = ifft( fft(y,NFFT) .* conj(fft(y,NFFT)) );

rt = roots( lpc(r, P) );
rt
rt = rt(imag(rt) > 0);
angz = atan2(imag(rt),real(rt));
[frqs, indices] = sort(angz.*(fs/(2*pi)), 'descend');

end

