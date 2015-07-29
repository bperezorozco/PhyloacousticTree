function [ dist ] = spectrogram_distance( filename1, filename2, samples )
%SPECTROGRAM_DISTANCE Summary of this function goes here
%   Detailed explanation goes here

x1 = audioread( filename1 );
x2 = audioread( filename2 );

min_length = min( length(x1), length(x2) );

if min_length < samples
    samples = min_length;
end

S1 = spectrogram( mean_normalise( x1(1:samples) ) );
S2 = spectrogram( mean_normalise( x2(1:samples) ) );

dist = mean( mean( log( abs( S1 - S2 ) ) ) );
end

