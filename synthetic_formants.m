y = [A*sin(2*pi*f*t) A3*sin(2*pi*f3*t) A2*sin(2*pi*f2*t)];
spectrogram(y, 100, 90, 256, fs, 'yaxis' );
F = formants(y', fs, 10);
figure
hold on
plot(F(:, 1))
plot(F(:, 2), 'r')
plot(F(:, 3), 'g')
plot(F(:, 5), 'k')
plot(F(:, 4), 'm')