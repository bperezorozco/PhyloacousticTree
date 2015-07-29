function [ Mel_spectra ] = mel_spectra( signal, fs, frame_length, overlap_length )
fft_size = 512;
fft_take = round( fft_size ) / 2 + 1;

end