function [ F ] = formants_from_mfcc_file( filename, n_formants, timeframe )
%FORMANTS_FILE Summary of this function goes here
%   Detailed explanation goes here
    %filename
    [x fs] = audioread( filename );
    F = formants_from_mfcc( x, fs, n_formants, timeframe );
end