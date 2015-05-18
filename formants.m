function [ F ] = formants( signal, fs, window_length, overlap, n_lpc, n_formants, Freqs, window, Emph )
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
%       integer, default 256 samples
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

if nargin < 2
    display( 'REQUIRED: @signal, @fs' );
    return;
end

if nargin < 9
    Emph = [1 0.63];
end
if nargin < 8
    window = 'hamming';
end
if nargin < 7
    Freqs = [ 90 400 ];
end
if nargin < 6
    n_formants = 4;
end
if nargin < 5
    n_lpc = 2 + round( fs / 1000 );%justify
end
if nargin < 4
    overlap = 40;
end
if nargin < 3
    window_length = 256;
end

shift = window_length - overlap;
len = length( signal ) - mod( length(signal), shift ) - shift;
f_threshold = 0;
frame = 1;

display( 'Working...' );

for begin = 1:shift:len
    x = signal( begin:(begin + window_length - 1) );
    
    if strcmp( window, 'hamming' )
        x = x .* hamming( length(x) );
    end
    
    if length( Emph ) == 2 && Emph(1) > 0
        x = filter(1, Emph, x);
    end
    
    a = lpc( x, n_lpc );
    rts = roots( a );
    rts = rts( imag(rts) >= f_threshold );
    angz = atan2( imag(rts), real(rts) );

    [frqs, indices] = sort( angz .* ( fs / ( 2 * pi ) ), 'descend' );
    bw = ( -1 / 2 ) * ( fs / ( 2 * pi ) ) * log( abs ( rts( indices ) ) );
    
    i_formants = 1;
    for i = 1:length(frqs)
        if ( frqs(i) > Freqs(1) && bw(i) < Freqs(2) )
            F(frame, i_formants) = frqs(i);
            i_formants = i_formants + 1;
        end
    end
    
    frame = frame + 1;
    
end

F = F(:, 1:n_formants);

display( 'Finished.' );
end



