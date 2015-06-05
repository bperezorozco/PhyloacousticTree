A = 1;
A2 = 3;
A3 = 5;
fs = 40000;
f = 2000;
f2 = 5000;
f3 = 7000;
T = fs;
t = (0:T-1)/fs;

y = A*sin(2*pi*f*t) + A2*sin(2*pi*f2*t) + A3*sin(2*pi*f3*t);
spectrogram(y, 100, 90, 256, fs, 'yaxis' );
[F B] = formants(y', fs, 10);
% figure
% hold on
% plot(F(:, 1))
% plot(F(:, 2), 'r')
% plot(F(:, 3), 'g')
% plot(F(:, 5), 'k')
% plot(F(:, 4), 'm')