function [ x ] = analyse_spectrograms( folders )
%ANALYSE_SPECTROGRAMS Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length( folders )
    foldername = strcat( './dataset/', folders{i} );
    files = dir( strcat( foldername, '/*.wav' ) );
    
    for j=1:length(files)
        filename = files(j).name;
        [x fs] = audioread( strcat(foldername, '/', filename) );
        spectrogram( mean_normalise(x), 256, 40, 256, fs, 'yaxis' );
        saveas( gcf, strcat('./plots/', filename, '.png') );
    end
end

end

