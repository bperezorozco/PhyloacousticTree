function [ output_args ] = vis_spectrogram( filename )
%VIS_SPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here
[x fs] = audioread( filename );
S = abs( spectrogram( mean_normalise( x ) ) );
figure;
imagesc(S);
end

