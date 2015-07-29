function [ F ] = formants_from_mfcc_file( filename, n_formants, timeframe )
%FORMANTS_FILE Summary of this function goes here
%   Detailed explanation goes here
    %filename
    [x fs] = audioread( filename );
    r = 2;
    F = formants_from_mfcc( resample( x, 1, r ), fs / r, n_formants, timeframe );
end