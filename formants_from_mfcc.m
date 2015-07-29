function [ F B ] = formants_from_mfcc( signal, fs, n_formants, window_length, overlap, n_lpc, Freqs, window, Emph )
%FORMANTS Estimates the formants of @signal every @window_length samples
%with @overlap samples. It returns the first @n_formants in the interval
%@Freqs. It can do preprocessing by applying window @window and
%pre-emphasis with coefficients @Emph
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REQUIRED ARGUMENTS:
%   @signal - vector with the signal to be analysed
%   @fs - sampling frequency

%DEFAULT PARAMETERS:
%   @window_length
%       integer, default 10 ms
%   @overlap 
%       integer, default 40 samples
%   @n_formants 
%       integer, default 2 formants
%   @Freqs 
%       vector of length 2 integers (in Hz), interval in which the formants
%           will be considered
%   @window 
%       default no window, otherwise 'Hamming'
%   @Emph 
%       vector of 2 decimals, coefficients for preemphasis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PAUSE = true;
DRAW_RESPONSE = 20;

if nargin < 3
    display( 'REQUIRED: @signal, @fs, @window_length' );
    return;
end

if nargin < 9
    Emph = [1 0.63];
end
if nargin < 8
    window = 'hamming';
end
if nargin < 7
    Freqs = [ 90 7000 400 ];
end
if nargin < 6
    n_lpc = n_formants * 2 + 2;
end
if nargin < 5
    overlap = 40;
end
if nargin < 4
    window_length = 10;
end
if nargin < 3
    n_formants = floor( fs / 2000 );
end


window_length = round( fs * window_length / 1000 );
overlap = round( window_length * 0.25 );
shift = window_length - overlap;
len = length( signal ) - mod( length(signal), shift ) - shift;
frames = len / shift;
F = zeros(frames, n_formants);
B = zeros(frames, n_formants);
threshold = 200;

frame = 1;
P = 2*n_formants + 2;

for begin = 1:shift:len
    y = signal( begin:(begin + window_length - 1) );
    if strcmp( window, 'hamming' )
        y = y .* hamming( length(y) );
    end
    
    if length( Emph ) == 2 && Emph(1) > 0
        y = filter(1, Emph, y);
    end
    
   % y = mean_normalise( y );
    NFFT = 2^nextpow2( length(y) );
    M = mel2freq( freq2mel( fs / 2 ) * linspace(0, 1, NFFT / 2) );
    f = fs / 2 * ( linspace(0, 1, NFFT / 2) );
    Y = fft(y, NFFT);
    r = ifft( Y .* conj( Y ) );
    
    mean_normalise( r );
    
    %a = lpc(r, P);
    a = levinson(r, P);
    if ismember( 1, isnan(a) )
        display('Found silent portion, skipping');
        continue;
    end
    
    rt = roots( a );
    rt = rt(imag(rt) >= 0);
    angz = atan2(imag(rt),real(rt));
    frqs = angz.*(fs/(2*pi));
    bw = -1/2*(fs/(2*pi))*log(abs(rt));
%     [bw, indices] = sort(bw, 'ascend');
%     frqs = frqs(indices);

    [frqs, indices] = sort(frqs, 'ascend');
    bw = bw(indices);
    frqs( frqs < Freqs(1) | frqs >= Freqs(2) |  bw > Freqs(3) ) = 0;
    [tmp1 tmp2] = filter_formants( frqs, threshold, n_formants, bw );
    F(frame, :) = tmp1;
    B(frame, :) = tmp2;
    
    if PAUSE && DRAW_RESPONSE == frame 
            [h,f]=freqz(1,a,512,fs);
            figure;
            plot(f,20*log10(abs(h)+eps));
            legend('LP Filter');
            xlabel('Frequency (Hz)');
            ylabel('Gain (dB)');
            DRAW_RESPONSE = DRAW_RESPONSE - 1;
    end
    
    frame = frame+1;
end

if PAUSE
    figure;
    NFFT = 2 ^ nextpow2( window_length );
    spectrogram(mean_normalise(signal), window_length, overlap, NFFT, fs, 'yaxis' );
    colormap bone;
    figure;
    plot(F(:, 1), 'b.');
    hold on;
    %plot(F(:, 2), 'r.');
    %plot(F(:, 3), 'g.');
    %figure;
    %plot(F(:, 1), F(:, 2), '.');
end

end



