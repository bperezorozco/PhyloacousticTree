function [ F ] = formants_file( filename )
%FORMANTS_FILE Summary of this function goes here
%   Detailed explanation goes here
    timeframe = 10;
    [x fs] = audioread( filename );
    F = formants( x, fs, timeframe );
end

